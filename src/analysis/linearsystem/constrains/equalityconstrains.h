
#ifndef SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_EQUALITYCONSTRAINS_H_
#define SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_EQUALITYCONSTRAINS_H_

#include "feti/feti.h"
#include "analysis/linearsystem/matrices/vector_distributed.h"

namespace espreso {

template <typename T>
struct EqualityConstrains {

	EqualityConstrains(FETI<T> &feti): feti(feti) {}

	void set(step::Step &step, const Vector_Distributed<Vector_Sparse, T> &dirichlet);
	void update(step::Step &step, const Vector_Distributed<Vector_Sparse, T> &dirichlet);

protected:
	FETI<T> &feti;
	std::vector<size_t> doffset;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_EQUALITYCONSTRAINS_H_ */
