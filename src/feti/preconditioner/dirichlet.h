
#ifndef SRC_FETI_PRECONDITIONER_DIRICHLET_H_
#define SRC_FETI_PRECONDITIONER_DIRICHLET_H_

#include "preconditioner.h"
#include "math/wrappers/math.spsolver.h"

namespace espreso {

template <typename T>
struct Dirichlet: public Preconditioner<T> {
    Dirichlet(FETI<T> &feti);
    virtual ~Dirichlet();

    void info() override;
    void update(const step::Step &step) override;

    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) override;

protected:
    void _manual(Matrix_CSR<T> &K, Matrix_Dense<T> &sc, std::vector<int> &permutation);
    void _print(const step::Step &step);

    using Preconditioner<T>::feti;
    std::vector<Vector_Dense<T> > Btx, KBtx;
    std::vector<Matrix_Dense<T> > sc;
    std::vector<std::vector<int> > indices;
};

}

#endif /* SRC_FETI_PRECONDITIONER_DIRICHLET_H_ */
