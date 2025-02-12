
#ifndef SRC_PHYSICS_KERNELS_BASEFUNCTIONS_VOLUME_HEXAHEDRON_H_
#define SRC_PHYSICS_KERNELS_BASEFUNCTIONS_VOLUME_HEXAHEDRON_H_

#include <vector>

namespace espreso {

struct Hexahedron {
    static int maxorder();
    static bool gpw(int order, std::vector<double> &r, std::vector<double> &s, std::vector<double> &t, std::vector<double> &w);
};

}

#endif /* SRC_PHYSICS_KERNELS_BASEFUNCTIONS_VOLUME_HEXAHEDRON_H_ */
