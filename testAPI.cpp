
#include "libespreso/feti4i.h"
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

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

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	int MPIrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

	std::stringstream path;
	path << argv[1] << "/" << MPIrank << "/";

	FETI4IMatrix K;
	FETI4IRHS rhs;

	FETI4ICreateStiffnessMatrixAndRHS(&K, &rhs, 0);

	std::vector<FETI4IInt> elements;
	readFile(elements, path.str() + "elements.txt");
	for (size_t e = 0; e < elements.size(); e++) {
		std::vector<FETI4IInt> indices;
		std::vector<FETI4IReal> Kvalues, RHSvalues;
		std::stringstream Ki, Kv, Rv;
		Ki << path.str() << "Ki" << elements[e] << ".txt";
		Kv << path.str() << "Ke" << elements[e] << ".txt";
		Rv << path.str() << "Rv" << elements[e] << ".txt";
		readFile(indices, Ki.str());
		readFile(Kvalues, Kv.str());
		readFile(RHSvalues, Rv.str());

		FETI4IAddElement(K, rhs, indices.size(), indices.data(), Kvalues.data(), RHSvalues.data());
	}

	std::vector<FETI4IInt> dirichlet_indices;
	std::vector<FETI4IReal> dirichlet_values;
	std::vector<FETI4IInt> l2g;
	std::vector<FETI4IInt> neighbours;
	readFile(dirichlet_indices, path.str() + "dirichlet_indices.txt");
	readFile(dirichlet_values, path.str() + "dirichlet_values.txt");
	readFile(l2g, path.str() + "l2g.txt");
	readFile(neighbours, path.str() + "neighbours.txt");

	FETI4IInstance instance;

	std::cout << "create instance\n";

	FETI4ICreateInstance(
			&instance,
			K,
			rhs,
			dirichlet_indices.size(),
			dirichlet_indices.data(),
			dirichlet_values.data(),
			l2g.data(),
			neighbours.size(),
			neighbours.data());

	std::vector<FETI4IReal> solution(l2g.size());

	std::cout << "solve\n";

	FETI4ISolve(instance, solution.size(), solution.data());



	MPI_Finalize();
}



