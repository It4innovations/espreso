
#ifndef SRC_FETI_PRECONDITIONER_DIRICHLET_IMPLICIT_H_
#define SRC_FETI_PRECONDITIONER_DIRICHLET_IMPLICIT_H_

#include "preconditioner.h"
#include "math/wrappers/math.spsolver.h"

namespace espreso {

template <typename T>
struct DirichletImplicit: public Preconditioner<T> {
    DirichletImplicit(FETI<T> &feti);
    virtual ~DirichletImplicit();

    void info() override;
    void set(const step::Step &step) override;
    void update(const step::Step &step) override;

    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) override;

protected:
    void _print(const step::Step &step);

    using Preconditioner<T>::feti;
    std::vector<Vector_Dense<T> > Btx, KBtx;
    std::vector<Matrix_Dense<T> > sc;
    std::vector<std::vector<int> > indices, permutation;

    std::vector<Matrix_CSR<T> > A11, A12, A22;
    std::vector<DirectSparseSolver<T, int> > solver;
};

}



#endif /* SRC_FETI_PRECONDITIONER_DIRICHLET_IMPLICIT_H_ */
