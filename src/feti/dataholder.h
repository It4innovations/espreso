
#ifndef SRC_ASSEMBLER_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_H_

#include "feti/generic/SparseMatrix.h"
#include "config/ecf/linearsolver/feti.h"
#include "math/data.decomposition.h"

#include <cstddef>
#include <vector>
#include <functional>

namespace espreso {

class SparseMatrix;

struct DataHolder {

	DataHolder();

	const DataDecomposition *decomposition;

	std::vector<SparseMatrix> origK, K, origKN1, origKN2, N1, N2, RegMat;
	std::vector<std::vector<double> > F;

	// matrices for Hybrid FETI constraints
	std::vector<SparseMatrix> B0;

	// matrices for FETI constraints
	std::vector<SparseMatrix> B1;
	std::vector<esint> B1Map;

	std::vector<std::vector<double> > B1c, LB, B1duplication;

	std::vector<std::vector<double> > primalSolution;
	std::vector<std::vector<double> > dualSolution;

	void assembleB0fromKernels();
};

}


#endif /* SRC_ASSEMBLER_INSTANCE_H_ */
