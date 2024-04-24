
#ifndef SRC_FETI_DUALOPERATOR_HYBRIDFETI_IMPLICIT_H_
#define SRC_FETI_DUALOPERATOR_HYBRIDFETI_IMPLICIT_H_

#include "dualoperator.h"
#include "math/wrappers/math.solver.h"
#include "math/wrappers/math.spsolver.h"

namespace espreso {

/*
 * K+: KxK : block diagonal
 * B1: LxK : from primal to dual
 * B0: CxK : from primal to cluster
 *
 * F0 G0t  u  =  d
 * G0  0   b  =  e
 *
 * F0 = B0 * K+ * B0t
 * G0 =     -Rt * B0t
 * d  =  B0 * K+ * B1 * x
 * e  =      -RT * B1 * x
 *
 * y = F * x = (B1 * K+ * B1t) * x + (B1 * K+ * B0) * u + (B1 * R * b)
 *
 */

template <typename T>
class HybridFETIImplicit: public DualOperator<T> {
public:
    HybridFETIImplicit(FETI<T> &feti);
    ~HybridFETIImplicit();

    void info();
    void set(const step::Step &step);
    void update(const step::Step &step);

    // y = (B1 * K+ * B1t) * x + (B1 * K+ * B0) * u + (B1 * R * b)
    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);

    // y = K+(f - Bt * x)
    void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y);

protected:
    void reduceInfo(DualOperatorInfo &sum, DualOperatorInfo &min, DualOperatorInfo &max);
    void printInfo(DualOperatorInfo &sum, DualOperatorInfo &min, DualOperatorInfo &max);

    void _applyK(const std::vector<Vector_Dense<T> > &x, std::vector<Vector_Dense<T> > &y);

    using DualOperator<T>::feti;
    using DualOperator<T>::d;

    std::vector<Matrix_CSR<T> > Kplus;
    std::vector<Vector_Dense<T> > Btx, KplusBtx;
    std::vector<DirectSparseSolver<T> > KSolver;

    std::vector<std::vector<int> > permutation;
    std::vector<Matrix_Dense<T> > dB0, dKB0, dF0;
    std::vector<int> G0offset;
    Matrix_CSR<T> F0, G0;
    DirectSparseSolver<T> F0Solver;
    DenseSolver<T> Splus;
    Vector_Dense<T> g, beta, mu;
};

}

#endif /* SRC_FETI_DUALOPERATOR_HYBRIDFETI_IMPLICIT_H_ */
