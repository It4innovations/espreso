
#ifndef SRC_FETI_PRECONDITIONER_PRECONDITIONER_H_
#define SRC_FETI_PRECONDITIONER_PRECONDITIONER_H_

#include "feti/feti.h"
#include "math/feti/vector_dual.h"

namespace espreso {

template <typename T>
class Preconditioner {
public:
	static Preconditioner<T>* set(FETI<T> *feti);

	Preconditioner(FETI<T> *feti): feti(feti) {}
	virtual ~Preconditioner() {}

	virtual void info() =0;
	virtual void update() =0;

	// y = S * x
	virtual void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) =0;

	FETI<T> *feti;
};

}

#endif /* SRC_FETI_PRECONDITIONER_PRECONDITIONER_H_ */
