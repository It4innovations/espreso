
#ifndef SRC_AXFETI_PRECONDITIONER_DIRICHLET_H_
#define SRC_AXFETI_PRECONDITIONER_DIRICHLET_H_

#include "preconditioner.h"

namespace espreso {

template <typename T>
class Dirichlet: public Preconditioner<T> {
public:
	Dirichlet(AX_FETI<T> *feti): Preconditioner<T>(feti) {}

	void info();
	void update();
};

}

#endif /* SRC_AXFETI_PRECONDITIONER_DIRICHLET_H_ */
