
#ifndef SRC_AXFETI_DUALOPERATOR_DUALOPERATOR_H_
#define SRC_AXFETI_DUALOPERATOR_DUALOPERATOR_H_

#include "axfeti/feti.h"

namespace espreso {

template <typename T>
class DualOperator {
public:
	static DualOperator<T>* set(AX_FETI<T> *feti);

	DualOperator(AX_FETI<T> *feti): feti(feti) {}
	virtual ~DualOperator() {}

	virtual void info() =0;
	virtual void update() =0;

	AX_FETI<T> *feti;

	Matrix_FETI<Matrix_CSR, T> Kplus;
};

}

#endif /* SRC_AXFETI_DUALOPERATOR_DUALOPERATOR_H_ */
