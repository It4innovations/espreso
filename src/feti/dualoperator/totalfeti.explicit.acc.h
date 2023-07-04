
#ifndef SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_ACC_H_
#define SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_ACC_H_

#include "totalfeti.explicit.h"

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
class TotalFETIExplicitAcc: public TotalFETIExplicit<T> {
public:
	TotalFETIExplicitAcc(FETI<T> &feti);
	~TotalFETIExplicitAcc();

	void info();
	void set(const step::Step &step);
	void update(const step::Step &step);

	// y = F * x
	void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);
	// y = K+(f - Bt * x)
	void toPrimal(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y);

protected:
	using DualOperator<T>::feti;
	using DualOperator<T>::d;
	using TotalFETIImplicit<T>::Kplus;
	using TotalFETIImplicit<T>::KSolver;

	std::vector<Matrix_CSC<T> > L, U;
	std::vector<Vector_Dense<int> > p;
};

}

#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_ACC_H_ */
