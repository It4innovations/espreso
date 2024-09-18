
#ifndef SRC_FETI_ITERATIVESOLVER_ORTHOCPG_H_
#define SRC_FETI_ITERATIVESOLVER_ORTHOCPG_H_

#include "iterativesolver.h"
#include "feti/common/matrix_dual_orthogonal.h"

namespace espreso {

template <typename T>
class OrthogonalizedCPG: public IterativeSolver<T> {
public:
    OrthogonalizedCPG(FETI<T> &feti);

    void info();

    void update(const step::Step &step)
    {
        IterativeSolver<T>::resize(l, r, w, x);
        IterativeSolver<T>::resize(pi, Fpi);
    }

    void solve(const step::Step &step, IterativeSolverInfo &info);

    using IterativeSolver<T>::feti;
    Vector_Dual<T> l, r, w, x;
    Matrix_Dual_Orthogonal<T> pi, Fpi;
    std::vector<T> wFp, pFp;
};

}


#endif /* SRC_FETI_ITERATIVESOLVER_ORTHOCPG_H_ */
