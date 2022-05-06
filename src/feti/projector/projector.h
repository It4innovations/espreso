
#ifndef SRC_FETI_PROJECTOR_PROJECTOR_H_
#define SRC_FETI_PROJECTOR_PROJECTOR_H_

#include "feti/feti.h"

#include "math/feti/vector_dual.h"
#include "math/feti/vector_kernel.h"

namespace espreso {

template <typename T>
class Projector {
public:
	static Projector<T>* set(FETI<T> *feti);

	Projector(FETI<T> *feti): feti(feti) {}
	virtual ~Projector() {}

	virtual void info() =0;
	virtual void update() =0;

	virtual void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) =0;
	virtual void applyGtInvGGt(const Vector_Kernel<T> &x, Vector_Dual<T> &y) =0;
	virtual void applyRInvGGtG(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y) =0;

	FETI<T> *feti;

	Vector_Kernel<T> e;
};

}

#endif /* SRC_FETI_PROJECTOR_PROJECTOR_H_ */
