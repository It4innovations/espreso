
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

// https://dl.acm.org/doi/pdf/10.1145/2929908.2929909

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
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y);

    // y = K+(f - Bt * x)
    void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y);
    void BtL(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y);

protected:
    using DualOperator<T>::feti;
    using DualOperator<T>::d;

    void _apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);
    void _applyK(std::vector<Vector_Dense<T> > &b, std::vector<Vector_Dense<T> > &x);

    std::vector<Matrix_CSR<T> > Kplus;
    std::vector<Vector_Dense<T> > Btx, KplusBtx;
    std::vector<DirectSparseSolver<T> > KSolver;

    void _computeB0();
    void _computeF0();
    void _computeG0();
    void _computeS0();

    void _compute_beta_mu(std::vector<Vector_Dense<T> > &b);

    std::vector<Matrix_CSR<T> > B0;
    std::vector<std::vector<int> > D2C, D2C0;
    std::vector<Matrix_Dense<T> > dKB0, origR1;
    std::vector<int> G0offset;
    Matrix_CSR<T> F0, G0;
    Matrix_Dense<T> S0;
    DirectSparseSolver<T> F0Solver;
    DenseSolver<T> Splus;
    Vector_Dense<T> g, beta, mu;
};

}

#endif /* SRC_FETI_DUALOPERATOR_HYBRIDFETI_IMPLICIT_H_ */
