
#ifndef SRC_FETI_DUALOPERATOR_DUALOPERATOR_H_
#define SRC_FETI_DUALOPERATOR_DUALOPERATOR_H_

#include "feti/feti.h"
#include "feti/common/vector_dual.h"

#include <vector>

namespace espreso {

struct DualOperatorInfo {
    size_t rows, nnzA, nnzL;
//     size_t memoryL;
    size_t dualA, surfaceA;
};

template <typename T>
class DualOperator {
public:
    static DualOperator<T>* set(FETI<T> &feti, const step::Step &step);

    DualOperator(FETI<T> &feti): feti(feti) {}
    virtual ~DualOperator() {}

    virtual void info() =0;
    virtual void set(const step::Step &step) =0;
    virtual void update(const step::Step &step) =0;

    // y = F * x
    virtual void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) =0;

    // y = K+(f - Bt * x)
    virtual void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y) =0;

    void estimateMaxEigenValue(double &lambda, int &iterations, double epsilon, int maxIterations);
    void estimateMaxProjectedEigenValue(double &lambda, int &iterations, double epsilon, int maxIterations, double rho = 1, double normPFP = 1);

    FETI<T> &feti;
    Vector_Dual<T> d;

protected:
    void getInitVector(Vector_Dual<T> &v);
};

}

#endif /* SRC_FETI_DUALOPERATOR_DUALOPERATOR_H_ */
