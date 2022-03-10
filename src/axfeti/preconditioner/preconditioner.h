
#ifndef SRC_AXFETI_PRECONDITIONER_PRECONDITIONER_H_
#define SRC_AXFETI_PRECONDITIONER_PRECONDITIONER_H_

#include "axfeti/feti.h"

namespace espreso {

template <typename T>
class Preconditioner {
public:
	static Preconditioner<T>* set(AX_FETI<T> *feti);

	Preconditioner(AX_FETI<T> *feti): feti(feti) {}
	virtual ~Preconditioner() {}

	virtual void info() =0;
	virtual void update() =0;

	AX_FETI<T> *feti;
};

}

#endif /* SRC_AXFETI_PRECONDITIONER_PRECONDITIONER_H_ */
