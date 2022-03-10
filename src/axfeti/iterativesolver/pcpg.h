
#ifndef SRC_AXFETI_ITERATIVESOLVER_PCPG_H_
#define SRC_AXFETI_ITERATIVESOLVER_PCPG_H_

#include "iterativesolver.h"

namespace espreso {

template <typename T>
class PCPG: public IterativeSolver<T> {
public:
	PCPG(AX_FETI<T> *feti);

	void info();
	void update();

//	void run(A, b, Proj, Prec)
//	{
//		w = Proj.leftapply(x);
//		z = A.apply(w);
//		zz = Prec.apply(z);
//		Proj.righ

//		evaluate(x);
//	}

};

}

#endif /* SRC_AXFETI_ITERATIVESOLVER_PCPG_H_ */
