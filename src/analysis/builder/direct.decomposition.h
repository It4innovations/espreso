
#ifndef SRC_ANALYSIS_BUILDER_DIRECT_DECOMPOSITION_H_
#define SRC_ANALYSIS_BUILDER_DIRECT_DECOMPOSITION_H_

#include <vector>

namespace espreso {

struct DirectDecomposition {
	esint begin, end, totalSize; // my DOFs
	std::vector<esint> halo; // halo indices
	std::vector<esint> neighDOF; // first DOF index per neighbor, the last is MY OFFSET
	std::vector<int> neighbors; // all neighboring process
};

}

#endif /* SRC_ANALYSIS_BUILDER_DIRECT_DECOMPOSITION_H_ */
