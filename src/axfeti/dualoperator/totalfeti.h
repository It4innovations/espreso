
#ifndef SRC_AXFETI_DUALOPERATOR_TOTALFETI_H_
#define SRC_AXFETI_DUALOPERATOR_TOTALFETI_H_

#include "dualoperator.h"

namespace espreso {

template <typename T>
class TotalFETI: public DualOperator<T> {
public:
	TotalFETI(AX_FETI<T> *feti);

	void info();
	void update();
};

}

#endif /* SRC_AXFETI_DUALOPERATOR_TOTALFETI_H_ */
