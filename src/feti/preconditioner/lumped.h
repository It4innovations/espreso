
#ifndef SRC_FETI_PRECONDITIONER_LUMPED_H_
#define SRC_FETI_PRECONDITIONER_LUMPED_H_

#include "preconditioner.h"

namespace espreso {

template <typename T>
class Lumped: public Preconditioner<T> {
public:
	Lumped(FETI<T> *feti): Preconditioner<T>(feti) {}
	~Lumped();

	void info() {}
	void update() {}

	void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);

	std::vector<Vector_Dense<T> > Btx, KBtx;
};

}

#endif /* SRC_FETI_PRECONDITIONER_LUMPED_H_ */
