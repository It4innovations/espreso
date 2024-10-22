
#ifndef SRC_ANALYSIS_PATTERN_DECOMPOSITION_DIRECT_H_
#define SRC_ANALYSIS_PATTERN_DECOMPOSITION_DIRECT_H_

#include <vector>

namespace espreso {

struct DecompositionDirect {
    esint begin, end, totalSize; // my DOFs
    std::vector<esint> halo; // halo indices
    std::vector<esint> neighDOF; // first DOF index per neighbor, the last is MY OFFSET
    std::vector<int> neighbors; // all neighboring process
};

}



#endif /* SRC_ANALYSIS_PATTERN_DECOMPOSITION_DIRECT_H_ */
