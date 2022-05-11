
#ifndef SRC_FETI_DUALOPERATOR_DUALOPERATOR_H_
#define SRC_FETI_DUALOPERATOR_DUALOPERATOR_H_

#include "feti/feti.h"
#include "math/feti/vector_dual.h"

#include <vector>

namespace espreso {

template <typename T>
class DualOperator {
public:
	static DualOperator<T>* set(FETI<T> *feti);

	DualOperator(FETI<T> *feti): feti(feti) {}
	virtual ~DualOperator() {}

	virtual void info() =0;
	virtual void update() =0;

	// y = F * x
	virtual void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) =0;

	// y = K+(f - Bt * x)
	virtual void toPrimal(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y) =0;

	FETI<T> *feti;

	std::vector<Matrix_CSR<T> > Kplus, KplusSurface;
	Vector_Dual<T> d;
};

}

#endif /* SRC_FETI_DUALOPERATOR_DUALOPERATOR_H_ */
