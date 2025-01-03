
#ifndef SRC_FETI_PRECONDITIONER_EMPTYPRECONDITIONER_H_
#define SRC_FETI_PRECONDITIONER_EMPTYPRECONDITIONER_H_

#include "preconditioner.h"
#include "math/math.h"

namespace espreso {

template <typename T>
class EmptyPreconditioner: public Preconditioner<T> {
public:
    EmptyPreconditioner(FETI<T> &feti): Preconditioner<T>(feti) {}

    void info() {}
    void update(const step::Step &step) {}

    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
    {
        math::copy(y, x);
    }

    virtual bool isset() const { return false; }
};

}

#endif /* SRC_FETI_PRECONDITIONER_EMPTYPRECONDITIONER_H_ */
