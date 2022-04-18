
#ifndef SRC_AXFETI_PROJECTOR_PROJECTOR_H_
#define SRC_AXFETI_PROJECTOR_PROJECTOR_H_

#include "axfeti/feti.h"

#include "math/feti/vector_dual.h"
#include "math/feti/vector_kernel.h"

namespace espreso {

template <typename T>
class Projector {
public:
	static Projector<T>* set(AX_FETI<T> *feti);

	Projector(AX_FETI<T> *feti): feti(feti) {}
	virtual ~Projector() {}

	virtual void info() =0;
	virtual void update() =0;

	virtual void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) =0;
	virtual void applyGtInvGGt(const Vector_Kernel<T> &x, Vector_Dual<T> &y) =0;
	virtual void applyInvGGtG(const Vector_Dual<T> &x, Vector_Kernel<T> &y) =0;

	AX_FETI<T> *feti;

	Vector_Kernel<T> e;
};

}

#endif /* SRC_AXFETI_PROJECTOR_PROJECTOR_H_ */
