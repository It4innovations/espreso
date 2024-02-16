
#ifndef SRC_FETI_ITERATIVESOLVER_CPG_H_
#define SRC_FETI_ITERATIVESOLVER_CPG_H_

#include "iterativesolver.h"

namespace espreso {

template <typename T>
class CPG: public IterativeSolver<T> {
public:
    CPG(FETI<T> &feti);

    void info();
    void solve(const step::Step &step, IterativeSolverInfo &info);

    using IterativeSolver<T>::feti;
    Vector_Dual<T> l, r, w, p;
    Vector_Dual<T> x, Fp;
};

}


#endif /* SRC_FETI_ITERATIVESOLVER_CPG_H_ */
