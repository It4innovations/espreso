
#ifndef SRC_ASSEMBLER_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_H_

#include <cstddef>
#include <vector>

namespace espreso {

class SparseMatrix;

struct NewInstance {

	NewInstance(size_t domains);

	size_t domains;
	std::vector<size_t> DOFs;

	std::vector<SparseMatrix> K, R1, R2, RegMat;
	std::vector<std::vector<double> > f;

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
	enum class CONSTRAINT {
		DIRICHLET,
		EQUALITY_CONSTRAINTS,
		INEQUALITY_CONSTRAINTS,
	};

	std::vector<size_t> block;


	std::vector<std::vector<double> > primalSolution;
	std::vector<std::vector<double> > dualSolution;
};

}



#endif /* SRC_ASSEMBLER_INSTANCE_H_ */
