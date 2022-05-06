
#ifndef SRC_FETI_PRECONDITIONER_EMPTYPRECONDITIONER_H_
#define SRC_FETI_PRECONDITIONER_EMPTYPRECONDITIONER_H_

#include "preconditioner.h"

namespace espreso {

template <typename T>
class EmptyPreconditioner: public Preconditioner<T> {
public:
	EmptyPreconditioner(FETI<T> *feti): Preconditioner<T>(feti) {}

	void info() {}
	void update() {}
};

}

#endif /* SRC_FETI_PRECONDITIONER_EMPTYPRECONDITIONER_H_ */
