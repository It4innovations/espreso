
#ifndef SRC_FETI_ITERATIVESOLVER_SMALSE_H_
#define SRC_FETI_ITERATIVESOLVER_SMALSE_H_

#include "iterativesolver.h"
#include "mprgp.h"

namespace espreso {

// Semi-Monotonic Augmented Lagrangian algorithm for Bound and Equality Constraints

// http://dx.doi.org/10.1007/s10107-011-0454-2

template <typename T>
class SMALBE: public IterativeSolver<T> {
public:
    SMALBE(FETI<T> &feti);

    void info();
    void solve(const step::Step &step, IterativeSolverInfo &info);

    using IterativeSolver<T>::feti;

    MPRGP<T> mprgp;
    Vector_Dual<T> b, b_, y, z, x_im, Fx_im, bCtmu, bCtmu_prev, Cx, Gx, gbCtmu;
    Vector_Kernel<T> mu;
};

}


#endif /* SRC_FETI_ITERATIVESOLVER_SMALSE_H_ */
