
#include "feti4i.h"

// we use ESPRESO assembler for compute stiffness matrices
#include "../src/input/meshgenerator/uniformmesh/cube/generator.h"
#include "../src/assembler/physics/linear/elasticity3d/assembler.h"

int main(int argc, char** argv)
{
	// Always initialize MPI before call ESPRESO!
	MPI_Init(&argc, &argv);

///////////////////////////////////////////// CREATE INSTANCE /////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	espreso::Mesh mesh; // for using ESPRESO to compute stiffness matrices we need to create mesh

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	espreso::input::CubeSettings settings(0, 1); // we assume only one MPI process
	settings.subdomainsInCluster[0] = settings.subdomainsInCluster[1] = settings.subdomainsInCluster[2] = 1;
	settings.elementsInSubdomain[0] = settings.elementsInSubdomain[1] = settings.elementsInSubdomain[2] = 6;
	settings.nodes["BOTTOM"] = espreso::Interval(0, settings.problemLength[0], 0, settings.problemLength[1], 0, 0);
	settings.properties["DIRICHLET"]["BOTTOM"] = "x: 0, y: 0, z: 0";

	// load mesh from cube generator
	espreso::input::CubeGenerator<espreso::input::Hexahedron8>::load(mesh, settings);

	// create assembler
	espreso::EqualityConstraints constraints(mesh);
	espreso::LinearElasticity3D assembler(mesh, constraints);
	assembler.prepareMeshStructures();
	assembler.assembleStiffnessMatrices();
	std::vector<espreso::Property> DOFs = espreso::LinearElasticity3D::pointDOFs;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<FETI4IInt> l2g(mesh.coordinates().clusterSize() * DOFs.size());
	for (size_t n = 0; n < mesh.coordinates().clusterSize(); n++) {
		for (size_t dof = 0; dof < DOFs.size(); dof++) {
			l2g[DOFs.size() * n + dof] = DOFs.size() * mesh.coordinates().globalIndex(n) + dof;
		}
	}

	std::vector<FETI4IInt> dirichletIndices;
	std::vector<FETI4IReal> dirichletValues;
	for (size_t n = 0; n < mesh.nodes().size(); n++) {
		for (size_t dof = 0; dof < DOFs.size(); dof++) {
			if (mesh.nodes()[n]->settings().isSet(DOFs[dof])) {
				dirichletIndices.push_back(DOFs.size() * n + dof);
				dirichletValues.push_back(mesh.nodes()[n]->settings(DOFs[dof]).back()->evaluate(n));
			}
		}
	}

	std::vector<FETI4IMPIInt> neighbours(mesh.neighbours().begin(), mesh.neighbours().end());
///////////////////////////////////////////// USE ESPRESO API /////////////////////////////////////////////////////////

	// At first create stiffness matrix
	FETI4IMatrix K;
	FETI4IInt indexBase = 0;
	FETI4ICreateStiffnessMatrix(&K, indexBase);

	// RHS: here, we assume that all nodes have the same number of DOFs, but it is not necessary
	std::vector<double> rhs(mesh.nodes().size() * DOFs.size());

	// Compose the matrix from elements matrices
	for (size_t e = 0; e < mesh.elements().size(); e++) {
		FETI4IInt dimension = static_cast<FETI4IInt>(mesh.elements()[e]->type());
		std::vector<FETI4IInt> nodes;
		std::vector<FETI4IInt> dofs;
		espreso::DenseMatrix Ke;
		std::vector<double> fe;

		assembler.assembleStiffnessMatrix(mesh.elements()[e], Ke, fe, dofs);
		for (size_t n = 0, i = 0; n < mesh.elements()[e]->nodes(); n++) {
			nodes.push_back(mesh.elements()[e]->node(n));
			for (size_t dof = 0; dof < DOFs.size(); dof++, i++) {
				rhs[dofs[i]] += fe[i];
			}
		}

		FETI4IAddElement(K, dimension, nodes.size(), nodes.data(), dofs.size(), dofs.data(), Ke.values());
	}

	FETI4IInt iopts[FETI4I_INTEGER_OPTIONS_SIZE];
	FETI4IReal ropts[FETI4I_REAL_OPTIONS_SIZE];

	FETI4ISetDefaultIntegerOptions(iopts);
	FETI4ISetDefaultRealOptions(ropts);

	// Configure ESPRESO solver
	iopts[FETI4I_SUBDOMAINS] = 4;
	iopts[FETI4I_PRECONDITIONER] = 1;
	iopts[FETI4I_VERBOSE_LEVEL] = 1;
	iopts[FETI4I_TESTING_LEVEL] = 0;
	iopts[FETI4I_MEASURE_LEVEL] = 0;

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

	espreso::output::VTK vtk(mesh, "results");
	vtk.storeGeometry();
	vtk.storeValues("api_result", DOFs.size(), solution, espreso::output::Store::ElementType::NODES);

	// Remove data
	FETI4IDestroy(K);
	FETI4IDestroy(instance);

	MPI_Finalize();
}



