
#ifndef SRC_AXFETI_DUALOPERATOR_TOTALFETI_H_
#define SRC_AXFETI_DUALOPERATOR_TOTALFETI_H_

#include "dualoperator.h"

namespace espreso {

/*
 * K+: KxK : block diagonal
 * B : LxK : from primal to dual
 *
 * y = F * x = B * K+ * Bt * x
 *
 * Btx = Bt * x          :: x        -> Btx     : L -> K (per domain)
 * KplusBtx = K+ * Btx   :: Btx      -> KplusBtx: K -> K (per domain)
 * y = B * KplusBtx      :: KplusBtx -> y       : K -> L (per domain except RBM and mortars)
 *
 */

template <typename T>
class TotalFETI: public DualOperator<T> {
public:
	TotalFETI(AX_FETI<T> *feti);
	~TotalFETI();

	void info();
	void update();

	// y = F * x
	void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);
	// y = K+(f - Bt * x)
	void toPrimal(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y);

	std::vector<Vector_Dense<T> > Btx, KplusBtx;
};

}

#endif /* SRC_AXFETI_DUALOPERATOR_TOTALFETI_H_ */
