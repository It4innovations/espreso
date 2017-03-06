
#include "feti4i.h"

// we use ESPRESO assembler for compute stiffness matrices
#include "../src/app/factory/factory.h"
#include "../src/assembler/instance/instance.h"
#include "../src/assembler/solution.h"
#include "../src/assembler/step.h"
#include "../src/basis/matrices/denseMatrix.h"
#include "../src/mesh/settings/property.h"
#include "../src/mesh/elements/element.h"
#include "../src/mesh/structures/mesh.h"
#include "../src/mesh/structures/elementtypes.h"
#include "../src/mesh/structures/coordinates.h"
#include "../src/assembler/old_physics/assembler.h"
#include "../src/configuration/globalconfiguration.h"
#include "../src/output/resultstore/vtkxmlascii.h"


int main(int argc, char** argv)
{
	// Always initialize MPI before call ESPRESO!
	MPI_Init(&argc, &argv);

///////////////////////////////////////////// CREATE INSTANCE /////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	espreso::GlobalConfiguration configuration(&argc, &argv);

	ESINFO(espreso::OVERVIEW) << "Run ESPRESO API test on " << espreso::environment->MPIsize << " process(es).";

	// use ESPRESO factory to allow run test on all examples/assemblers
	espreso::Factory factory(configuration);
	factory.mesh->partitiate(1); // we need only one sub-domain per cluster
	factory.instance->init(); // it perform full initialization; it is not effective but sufficient to test API

	std::vector<espreso::Property> DOFs = factory.instance->physics().pointDOFs;
	size_t DOFsSize = factory.instance->physics().f[0].size();

	std::vector<FETI4IInt> l2g(DOFsSize);

	// TODO: generalize l2g and dirichlet
	for (size_t n = 0; n < factory.mesh->coordinates().clusterSize(); n++) {
		for (size_t dof = 0; dof < DOFs.size(); dof++) {
			l2g[DOFs.size() * n + dof] = DOFs.size() * factory.mesh->coordinates().globalIndex(n) + dof;
		}
	}

	std::vector<FETI4IInt> dirichletIndices;
	std::vector<FETI4IReal> dirichletValues;
	for (size_t n = 0; n < factory.mesh->nodes().size(); n++) {
		for (size_t dof = 0; dof < DOFs.size(); dof++) {
			if (factory.mesh->nodes()[n]->hasProperty(DOFs[dof], 0)) {
				dirichletIndices.push_back(DOFs.size() * n + dof);
				dirichletValues.push_back(factory.mesh->nodes()[n]->getProperty(DOFs[dof], 0, 0, 0));
			}
		}
	}
/////////////////////////////////////////////// USE ESPRESO API /////////////////////////////////////////////////////////

	// At first create stiffness matrix
	FETI4IMatrix K;
	FETI4IInt indexBase = 0;
	FETI4ICreateStiffnessMatrix(&K, (int)factory.instance->physics().mtype, indexBase);

	std::vector<FETI4IReal> rhs(factory.instance->physics().f[0]);
	std::vector<FETI4IMPIInt> neighbours(factory.mesh->neighbours());

	// Compose the matrix from elements matrices
	for (size_t e = 0; e < factory.mesh->elements().size(); e++) {
		FETI4IInt dimension = static_cast<FETI4IInt>(factory.mesh->elements()[e]->type());
		std::vector<FETI4IInt> nodes;
		std::vector<FETI4IInt> dofs;
		espreso::DenseMatrix Ke;
		std::vector<double> fe;

		for (size_t n = 0; n < factory.mesh->elements()[e]->nodes(); n++) {
			nodes.push_back(factory.mesh->elements()[e]->node(n));
		}

		factory.instance->physics().assembleStiffnessMatrix(factory.mesh->elements()[e], Ke, fe, dofs);
		FETI4IAddElement(K, dimension, nodes.size(), nodes.data(), dofs.size(), dofs.data(), Ke.values());
	}

	FETI4IInt iopts[FETI4I_INTEGER_OPTIONS_SIZE];
	FETI4IReal ropts[FETI4I_REAL_OPTIONS_SIZE];

	// FETI4ISetDefaultIntegerOptions(iopts);
	// FETI4ISetDefaultRealOptions(ropts);

	// Configure ESPRESO solver
	iopts[FETI4I_FETI_METHOD] = static_cast<int>(configuration.linear_elasticity_2D.espreso.method);
	iopts[FETI4I_CGSOLVER] = static_cast<int>(configuration.linear_elasticity_3D.espreso.solver);
	iopts[FETI4I_SUBDOMAINS] = configuration.esdata.domains;
	iopts[FETI4I_ITERATIONS] = configuration.linear_elasticity_3D.espreso.iterations;
	iopts[FETI4I_PRECONDITIONER] = static_cast<int>(configuration.linear_elasticity_3D.espreso.preconditioner);
	iopts[FETI4I_VERBOSE_LEVEL] = configuration.env.verbose_level;
	iopts[FETI4I_TESTING_LEVEL] = configuration.env.testing_level;
	iopts[FETI4I_MEASURE_LEVEL] = configuration.env.measure_level;
	iopts[FETI4I_PRINT_MATRICES] = configuration.env.print_matrices;
	ropts[FETI4I_EPSILON] = configuration.linear_elasticity_3D.espreso.epsilon;


	// Create instance of a problem
	FETI4IInstance instance;
	FETI4ICreateInstance(
			&instance,
			K,
			rhs.size(),
			rhs.data(),
			l2g.data(),
			neighbours.size(),
			neighbours.data(),
			dirichletIndices.size(),
			dirichletIndices.data(),
			dirichletValues.data(),
			iopts,
			ropts);

	// Prepare memory for save solution
	std::vector<std::vector<FETI4IReal> > solution(1, std::vector<FETI4IReal>(rhs.size()));

	// Solve the system
	FETI4ISolve(instance, solution[0].size(), solution[0].data());

	// Process solution

	espreso::output::VTKXMLASCII vtk(configuration.output, factory.mesh, "results");
	espreso::Step step;
	std::vector<espreso::Solution*> sol;
	sol.push_back(new espreso::Solution("solution", espreso::ElementType::NODES, factory.instance->physics().pointDOFs, solution));
	vtk.storeSolution(step, sol);

	// Remove data
	FETI4IDestroy(K);
	FETI4IDestroy(instance);

	MPI_Finalize();
}



