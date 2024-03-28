
#ifndef SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_ELASTICITY_H_
#define SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_ELASTICITY_H_

#include "feti/feti.h"

namespace espreso {

template <typename T>
struct RegularizationElasticity {

    static void set(FETI<T> &feti);
    static void update(FETI<T> &feti);

    static void getFixPoints(std::vector<esint> &fixPoints, int domain, bool onSurface = false);
    static void getCorners(std::vector<esint> &fixPoints, int domain);

protected:
    static void set2D(FETI<T> &feti, esint domain);
    static void set3D(FETI<T> &feti, esint domain, bool onSurface);

    static std::vector<Matrix_Dense<T> > NtNNtN;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_ELASTICITY_H_ */
