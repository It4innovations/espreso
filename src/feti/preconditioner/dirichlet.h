
#ifndef SRC_FETI_PRECONDITIONER_DIRICHLET_H_
#define SRC_FETI_PRECONDITIONER_DIRICHLET_H_

#include "preconditioner.h"
#include "basis/containers/allocators.h"

namespace espreso {

template <typename T>
class Dirichlet: public Preconditioner<T> {
public:
	Dirichlet(FETI<T> *feti): Preconditioner<T>(feti) {}
	~Dirichlet();

	void info();
	void update();

	void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);

	std::vector<Vector_Dense<T> > Btx, KBtx;

	std::vector<esint, initless_allocator<esint> > surface;
};

}

#endif /* SRC_FETI_PRECONDITIONER_DIRICHLET_H_ */
