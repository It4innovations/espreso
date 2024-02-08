
#ifndef SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_CONSTRAINS_H_
#define SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_CONSTRAINS_H_

#include "analysis/math/vector_distributed.h"
#include "feti/feti.h"

namespace espreso {

template <typename T>
struct Constrains {

	static void set(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet);
	static void update(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet);
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_CONSTRAINS_H_ */
