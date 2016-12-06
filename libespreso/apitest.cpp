
#include "feti4i.h"

// we use ESPRESO assembler for compute stiffness matrices
#include "../src/app/factory/factory.h"
#include "../src/config/description.h"

int main(int argc, char** argv)
{
	// Always initialize MPI before call ESPRESO!
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

///////////////////////////////////////////// CREATE INSTANCE /////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	espreso::ParametersReader::fromArguments(&argc, &argv);

	ESINFO(espreso::OVERVIEW) << "Run ESPRESO API test on " << espreso::config::env::MPIsize << " process(es).";
	espreso::ParametersReader::printParameters(espreso::config::parameters, espreso::config::info::VERBOSE_LEVEL);

	// use ESPRESO factory to allow run test on all examples/assemblers
	espreso::Factory factory(espreso::configuration);
	factory.mesh.partitiate(1); // we need only one sub-domain per cluster
	factory.instance->init(); // it perform full initialization; it is not effective but sufficient to test API

	std::vector<espreso::Property> DOFs = factory.instance->physics().pointDOFs;
	size_t DOFsSize = factory.instance->physics().f[0].size();

	std::vector<FETI4IInt> l2g(DOFsSize);

	// TODO: generalize l2g and dirichlet
	for (size_t n = 0; n < factory.mesh.coordinates().clusterSize(); n++) {
		for (size_t dof = 0; dof < DOFs.size(); dof++) {
			l2g[DOFs.size() * n + dof] = DOFs.size() * factory.mesh.coordinates().globalIndex(n) + dof;
		}
	}

	std::vector<FETI4IInt> dirichletIndices;
	std::vector<FETI4IReal> dirichletValues;
	for (size_t n = 0; n < factory.mesh.nodes().size(); n++) {
		for (size_t dof = 0; dof < DOFs.size(); dof++) {
			if (factory.mesh.nodes()[n]->settings().isSet(DOFs[dof])) {
				dirichletIndices.push_back(DOFs.size() * n + dof);
				dirichletValues.push_back(factory.mesh.nodes()[n]->settings(DOFs[dof]).back()->evaluate(n));
			}
		}
	}
/////////////////////////////////////////////// USE ESPRESO API /////////////////////////////////////////////////////////

	// At first create stiffness matrix
	FETI4IMatrix K;
	FETI4IInt indexBase = 0;
	FETI4ICreateStiffnessMatrix(&K, (int)factory.instance->physics().mtype, indexBase);

	std::vector<FETI4IReal> rhs(factory.instance->physics().f[0]);
	std::vector<FETI4IMPIInt> neighbours(factory.mesh.neighbours());

	// Compose the matrix from elements matrices
	for (size_t e = 0; e < factory.mesh.elements().size(); e++) {
		FETI4IInt dimension = static_cast<FETI4IInt>(factory.mesh.elements()[e]->type());
		std::vector<FETI4IInt> nodes;
		std::vector<FETI4IInt> dofs;
		espreso::DenseMatrix Ke;
		std::vector<double> fe;

		for (size_t n = 0; n < factory.mesh.elements()[e]->nodes(); n++) {
			nodes.push_back(factory.mesh.elements()[e]->node(n));
		}

		factory.instance->physics().assembleStiffnessMatrix(factory.mesh.elements()[e], Ke, fe, dofs);
		FETI4IAddElement(K, dimension, nodes.size(), nodes.data(), dofs.size(), dofs.data(), Ke.values());
	}

	FETI4IInt iopts[FETI4I_INTEGER_OPTIONS_SIZE];
	FETI4IReal ropts[FETI4I_REAL_OPTIONS_SIZE];

	FETI4ISetDefaultIntegerOptions(iopts);
	FETI4ISetDefaultRealOptions(ropts);

	// Configure ESPRESO solver
	iopts[FETI4I_FETI_METHOD] = 0;
	iopts[FETI4I_SUBDOMAINS] = 4;
	iopts[FETI4I_ITERATIONS] = 100;
	iopts[FETI4I_PRECONDITIONER] = 1;
	iopts[FETI4I_VERBOSE_LEVEL] = 1;
	iopts[FETI4I_TESTING_LEVEL] = 0;
	iopts[FETI4I_MEASURE_LEVEL] = 0;
	iopts[FETI4I_PRINT_MATRICES] = 1;


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

	espreso::store::VTK vtk(factory.mesh, "results");
	vtk.storeGeometry();
	vtk.storeValues("api_result", DOFs.size(), solution, espreso::store::Store::ElementType::NODES);

	// Remove data
	FETI4IDestroy(K);
	FETI4IDestroy(instance);

	MPI_Finalize();
}



