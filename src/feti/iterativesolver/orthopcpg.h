
#ifndef SRC_FETI_ITERATIVESOLVER_ORTHOPCPG_H_
#define SRC_FETI_ITERATIVESOLVER_ORTHOPCPG_H_

#include "iterativesolver.h"
#include "feti/common/matrix_dual_orthogonal.h"

namespace espreso {

template <typename T>
class OrthogonalizedPCPG: public IterativeSolver<T> {
public:
    OrthogonalizedPCPG(FETI<T> &feti);

    void info() override;

    void update(const step::Step &step) override
    {
        IterativeSolver<T>::resize(l, r, w, y, z, x);
        IterativeSolver<T>::resize(pi, Fpi);
    }

    void solve(const step::Step &step, IterativeSolverInfo &info) override;

    using IterativeSolver<T>::feti;
    Vector_Dual<T> l, r, w, y, z, x;
    Matrix_Dual_Orthogonal<T> pi, Fpi;
    std::vector<T> yFp, pFp;
};

}

#endif /* SRC_FETI_ITERATIVESOLVER_ORTHOPCPG_H_ */
