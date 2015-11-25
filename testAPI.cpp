
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

static void readBinary(std::vector<double> &vector, std::string fileName) {
	std::ifstream file(fileName.c_str(), std::fstream::binary);
	if (file.is_open()) {
		for (size_t i = 0; i < vector.size(); i++) {
			double value;
			file.read(reinterpret_cast<char *>(&value), sizeof(double));
			vector[i] = value;
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

	FETI4ICreateStiffnessMatrix(&K, 0);

	std::vector<FETI4IInt> elements;
	readFile(elements, path.str() + "elements.txt");
	for (size_t e = 0; e < elements.size(); e++) {
		std::vector<FETI4IInt> indices;
		std::vector<FETI4IReal> Kvalues;
		std::stringstream Ki, Kv;
		Ki << path.str() << "Ki" << elements[e] << ".txt";
		Kv << path.str() << "Ke" << elements[e] << ".bin";
		readFile(indices, Ki.str());
		Kvalues.resize(indices.size() * indices.size());
		readBinary(Kvalues, Kv.str());

		FETI4IAddElement(K, indices.size(), indices.data(), Kvalues.data());
	}

	std::vector<FETI4IReal> rhs;
	std::vector<FETI4IInt> dirichlet_indices;
	std::vector<FETI4IReal> dirichlet_values;
	std::vector<FETI4IInt> l2g;
	std::vector<FETI4IMPIInt> neighbours;
	readFile(rhs, path.str() + "rhs.txt");
	readFile(dirichlet_indices, path.str() + "dirichlet_indices.txt");
	readFile(dirichlet_values, path.str() + "dirichlet_values.txt");
	readFile(l2g, path.str() + "l2g.txt");
	readFile(neighbours, path.str() + "neighbours.txt");

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

	std::vector<FETI4IReal> solution(rhs.size());

	FETI4ISolve(instance, solution.size(), solution.data());


	MPI_Finalize();
}



