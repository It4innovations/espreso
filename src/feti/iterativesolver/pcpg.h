
#ifndef SRC_FETI_ITERATIVESOLVER_PCPG_H_
#define SRC_FETI_ITERATIVESOLVER_PCPG_H_

#include "iterativesolver.h"

namespace espreso {

template <typename T>
class PCPG: public IterativeSolver<T> {
public:
    PCPG(FETI<T> &feti);

    void info();

    void update(const step::Step &step)
    {
        IterativeSolver<T>::resize(l, r, w, y, z, p);
        IterativeSolver<T>::resize(x, Fp);
    }

    void solve(const step::Step &step, IterativeSolverInfo &info);

    using IterativeSolver<T>::feti;
    Vector_Dual<T> l, r, w, y, z, p;
    Vector_Dual<T> x, Fp;
};

}

#endif /* SRC_FETI_ITERATIVESOLVER_PCPG_H_ */
