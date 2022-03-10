
#ifndef SRC_AXFETI_ITERATIVESOLVER_ITERATIVESOLVER_H_
#define SRC_AXFETI_ITERATIVESOLVER_ITERATIVESOLVER_H_

#include "axfeti/feti.h"

namespace espreso {

template <typename T>
class IterativeSolver {
public:
	static IterativeSolver<T>* set(AX_FETI<T> *feti);

	IterativeSolver(AX_FETI<T> *feti): feti(feti) {}
	virtual ~IterativeSolver() {}

	virtual void info() =0;
	virtual void update() =0;

	AX_FETI<T> *feti;
};

}

#endif /* SRC_AXFETI_ITERATIVESOLVER_ITERATIVESOLVER_H_ */
