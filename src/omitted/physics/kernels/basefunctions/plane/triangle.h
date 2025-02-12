
#ifndef SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_TRIANGLE_H_
#define SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_TRIANGLE_H_

#include <vector>

namespace espreso {

struct Triangle {
    static int maxorder();
    static bool gpw(int order, std::vector<double> &r, std::vector<double> &s, std::vector<double> &w);
};

}

#endif /* SRC_PHYSICS_KERNELS_BASEFUNCTIONS_PLANE_TRIANGLE_H_ */
