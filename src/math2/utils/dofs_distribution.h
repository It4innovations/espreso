
#ifndef SRC_MATH2_UTILS_DOFS_DISTRIBUTION_H_
#define SRC_MATH2_UTILS_DOFS_DISTRIBUTION_H_

#include <vector>

namespace espreso {

struct DOFsDistribution {
	esint begin, end, totalSize; // my DOFs
	std::vector<esint> halo; // halo indices
	std::vector<esint> neighDOF; // first DOF index per neighbor
	std::vector<int> neighbors; // all neighboring process
};

}

#endif /* SRC_MATH2_UTILS_DOFS_DISTRIBUTION_H_ */
