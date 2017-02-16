
#ifndef SRC_ASSEMBLER_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_H_

#include <cstddef>
#include <vector>
#include <fstream>

namespace espreso {

class SparseMatrix;
class Solution;

enum Matrices : int {
	K = 1,
	M = 2,
	R = 4,
	f = 8
};

struct Instance {

	Instance(size_t domains);
	~Instance();

	size_t domains;
	std::vector<size_t> DOFs;

	std::vector<SparseMatrix> K, M, N1, N2, RegMat;
	std::vector<std::vector<double> > R, f;

	// matrices for Hybrid FETI constraints
	std::vector<SparseMatrix> B0;
	std::vector<std::vector<esglobal> > B0subdomainsMap; // TODO: not needed

	// matrices for FETI constraints
	std::vector<SparseMatrix> B1;
	std::vector<std::vector<esglobal> > B1subdomainsMap; // TODO: not needed
	std::vector<std::vector<esglobal> > B1clustersMap; // TODO: get it directly

	std::vector<std::vector<double> > B1c;
	std::vector<std::vector<double> > LB;
	std::vector<std::vector<double> > B1duplicity;

	std::vector<SparseMatrix> inequality;
	std::vector<std::vector<double> > inequalityC;

	// blocks types of B1
	enum CONSTRAINT {
		DIRICHLET,
		EQUALITY_CONSTRAINTS,
		INEQUALITY_CONSTRAINTS,
	};

	std::vector<size_t> block;


	std::vector<std::vector<double> > primalSolution;
	std::vector<std::vector<double> > dualSolution;

	std::vector<Solution*> solutions;
};

inline Matrices operator|(Matrices m1, Matrices m2)
{
	return static_cast<Matrices>(static_cast<int>(m1) | static_cast<int>(m2));
}

inline Matrices operator|=(Matrices m1, Matrices m2)
{
	return static_cast<Matrices>(static_cast<int>(m1) | static_cast<int>(m2));
}

inline Matrices operator&(Matrices m1, Matrices m2)
{
	return static_cast<Matrices>(static_cast<int>(m1) & static_cast<int>(m2));
}

inline Matrices operator&=(Matrices m1, Matrices m2)
{
	return static_cast<Matrices>(static_cast<int>(m1) & static_cast<int>(m2));
}

}


#endif /* SRC_ASSEMBLER_INSTANCE_H_ */
