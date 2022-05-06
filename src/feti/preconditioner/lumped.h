
#ifndef SRC_FETI_PRECONDITIONER_LUMPED_H_
#define SRC_FETI_PRECONDITIONER_LUMPED_H_

#include "preconditioner.h"

namespace espreso {

template <typename T>
class Lumped: public Preconditioner<T> {
public:
	Lumped(FETI<T> *feti): Preconditioner<T>(feti) {}

	void info() {}
	void update() {}
};

}

#endif /* SRC_FETI_PRECONDITIONER_LUMPED_H_ */
