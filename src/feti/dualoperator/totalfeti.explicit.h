
#ifndef SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_H_
#define SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_H_

#include "dualoperator.h"
#include "math/wrappers/math.spsolver.h"

namespace espreso {

/*
 * K+: KxK : block diagonal
 * B : LxK : from primal to dual
 *
 * y = F * x = (B * K+ * Bt) * x
 *
 * Btx = Bt * x          :: x        -> Btx     : L -> K (per domain)
 * KplusBtx = K+ * Btx   :: Btx      -> KplusBtx: K -> K (per domain)
 * y = B * KplusBtx      :: KplusBtx -> y       : K -> L (per domain except RBM and mortars)
 *
 */

template <typename T>
class TotalFETIExplicit: public DualOperator<T> {
public:
    TotalFETIExplicit(FETI<T> &feti);
    virtual ~TotalFETIExplicit();

    void info() override;
    void set(const step::Step &step) override;
    void update(const step::Step &step) override;

    // y = F * x
    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) override;
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y) override;

    // y = K+(f - Bt * x)
    void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y) override;

protected:
    void _apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);

    void print(const step::Step &step);

    using DualOperator<T>::feti;
    using DualOperator<T>::d;

    std::vector<Matrix_CSR<T> > Kplus;
    std::vector<Vector_Dense<T> > Btx, KplusBtx;
    std::vector<DirectSparseSolver<T> > KSolver;

    std::vector<Matrix_Dense<T> > F;
    std::vector<Vector_Dense<T> > in, out;
};

}

#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_H_ */
