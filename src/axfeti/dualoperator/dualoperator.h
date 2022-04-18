
#ifndef SRC_AXFETI_DUALOPERATOR_DUALOPERATOR_H_
#define SRC_AXFETI_DUALOPERATOR_DUALOPERATOR_H_

#include "axfeti/feti.h"
#include "math/feti/vector_dual.h"

#include <vector>

namespace espreso {

template <typename T>
class DualOperator {
public:
	static DualOperator<T>* set(AX_FETI<T> *feti);

	DualOperator(AX_FETI<T> *feti): feti(feti) {}
	virtual ~DualOperator() {}

	virtual void info() =0;
	virtual void update() =0;

	// y = F * x
	virtual void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) =0;

	AX_FETI<T> *feti;

	std::vector<Matrix_CSR<T> > Kplus;
	Vector_Dual<T> d;
};

}

#endif /* SRC_AXFETI_DUALOPERATOR_DUALOPERATOR_H_ */
