
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

#include "feti4i.h"
#include "../src/config/esconfig.h"

template <typename Ttype>
static void readFile(typename std::vector<Ttype> &vector, std::string fileName) {
	vector.clear();
	std::ifstream file(fileName.c_str());
	if (file.is_open()) {
		Ttype value;
		while (file >> value) {
			vector.push_back(value);
		}
	} else {
		std::cerr << "Cannot read file " << fileName << "\n";
		exit(EXIT_FAILURE);
	}
}

static void readBinary(std::vector<FETI4IReal> &vector, std::string fileName) {
	std::ifstream file(fileName.c_str(), std::fstream::binary);
	if (file.is_open()) {
		for (size_t i = 0; i < vector.size(); i++) {
			FETI4IReal value;
			file.read(reinterpret_cast<char *>(&value), sizeof(double));
			vector[i] = value;
		}
	} else {
		std::cerr << "Cannot read file " << fileName << "\n";
		exit(EXIT_FAILURE);
	}
}


static void loadStructures(
		const std::string                    &path,

		std::vector<std::vector<FETI4IInt> > &K_indices,
		std::vector<std::vector<FETI4IReal> > &K_values,

		std::vector<FETI4IReal>              &rhs,

		std::vector<FETI4IInt>               &dirichlet_indices,
		std::vector<FETI4IReal>              &dirichlet_values,

		std::vector<FETI4IInt>               &l2g,

		std::vector<FETI4IMPIInt>            &neighbours
		)
{
	int MPIrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

	std::stringstream root;
	root << path << "/" << MPIrank << "/";

	std::vector<FETI4IInt> elements;
	readFile(elements, root.str() + "elements.txt");
	K_indices.resize(elements.size());
	K_values.resize(elements.size());
	for (size_t e = 0; e < elements.size(); e++) {
		std::stringstream Ki, Kv;
		Ki << root.str() << "Ki" << elements[e] << ".txt";
		Kv << root.str() << "Ke" << elements[e] << ".bin";
		readFile(K_indices[e], Ki.str());
		K_values[e].resize(K_indices[e].size() * K_indices[e].size());
		readBinary(K_values[e], Kv.str());
	}

	readFile(rhs,               root.str() + "rhs.txt");
	readFile(dirichlet_indices, root.str() + "dirichlet_indices.txt");
	readFile(dirichlet_values,  root.str() + "dirichlet_values.txt");
	readFile(l2g,               root.str() + "l2g.txt");
	readFile(neighbours,        root.str() + "neighbours.txt");

}


int main(int argc, char** argv)
{
	// Always initialize MPI before call ESPRESO!
	MPI_Init(&argc, &argv);

	// The following data are considered to be computed by a library
	std::vector<std::vector<FETI4IInt> >  K_indices;         // vector of stiffness matrices indices
	std::vector<std::vector<FETI4IReal> > K_values;          // vector of stiffness matrices values
	std::vector<FETI4IReal>               rhs;               // right hand side
	std::vector<FETI4IInt>                dirichlet_indices; // DOFs with dirichlet
	std::vector<FETI4IReal>               dirichlet_values;  // dirichlet values
	std::vector<FETI4IInt>                l2g;               // local to global numbering
	std::vector<FETI4IMPIInt>             neighbours;        // vector of neighbours

	// We load data stored in ESPRESO example directory
	loadStructures("../examples/api/cube", K_indices, K_values, rhs, dirichlet_indices, dirichlet_values, l2g, neighbours);


/*------------------------------------------------------------------------------
           When we have data, use ESPRESO API to compute results
------------------------------------------------------------------------------*/

	// At first create stiffness matrix
	FETI4IMatrix K;
	FETI4ICreateStiffnessMatrix(&K, 0);

	// Compose the matrix from elements matrices
	for (size_t i = 0; i < K_indices.size(); i++) {
		FETI4IAddElement(K, K_indices[i].size(), K_indices[i].data(), K_values[i].data());
	}

	// Configure ESPRESO solver
	espreso::config::mesh::SUBDOMAINS = 8;
	espreso::config::solver::PRECONDITIONER = espreso::config::solver::PRECONDITIONERalternative::DIRICHLET;

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
			dirichlet_indices.size(),
			dirichlet_indices.data(),
			dirichlet_values.data());

	// Prepare memory for save solution
	std::vector<FETI4IReal> solution(rhs.size());

	// Solve the system
	FETI4ISolve(instance, solution.size(), solution.data());

	// Process solution

	MPI_Finalize();
}






