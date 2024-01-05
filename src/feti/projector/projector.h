
#ifndef SRC_FETI_PROJECTOR_PROJECTOR_H_
#define SRC_FETI_PROJECTOR_PROJECTOR_H_

#include "feti/feti.h"
#include "feti/common/vector_dual.h"
#include "feti/common/vector_kernel.h"

namespace espreso {

template <typename T>
struct Projector {
	static Projector<T>* set(FETI<T> &feti, const step::Step &step);

	Projector(FETI<T> &feti): feti(feti) {}
	virtual ~Projector() {}

	virtual void info() =0;
	virtual void update(const step::Step &step) =0;

	virtual void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) =0;
	virtual void applyGtInvGGt(const Vector_Kernel<T> &x, Vector_Dual<T> &y) =0;
	virtual void applyRInvGGtG(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y) =0;

	Vector_Kernel<T> e;

protected:
	FETI<T> &feti;
};

}

#endif /* SRC_FETI_PROJECTOR_PROJECTOR_H_ */
