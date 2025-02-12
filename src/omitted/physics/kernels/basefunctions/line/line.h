
#ifndef SRC_PHYSICS_KERNELS_BASEFUNCTIONS_LINE_LINE_H_
#define SRC_PHYSICS_KERNELS_BASEFUNCTIONS_LINE_LINE_H_

#include <vector>

namespace espreso {

struct Line {
    static int maxorder();
    static bool gpw(int order, std::vector<double> &r, std::vector<double> &w);
};

}

#endif /* SRC_PHYSICS_KERNELS_BASEFUNCTIONS_LINE_LINE_H_ */
