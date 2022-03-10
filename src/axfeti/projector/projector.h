
#ifndef SRC_AXFETI_PROJECTOR_PROJECTOR_H_
#define SRC_AXFETI_PROJECTOR_PROJECTOR_H_

#include "axfeti/feti.h"

namespace espreso {

template <typename T>
class Projector {
public:
	static Projector<T>* set(AX_FETI<T> *feti);

	Projector(AX_FETI<T> *feti): feti(feti) {}
	virtual ~Projector() {}

	virtual void info() =0;
	virtual void update() =0;

	AX_FETI<T> *feti;
};

}

#endif /* SRC_AXFETI_PROJECTOR_PROJECTOR_H_ */
