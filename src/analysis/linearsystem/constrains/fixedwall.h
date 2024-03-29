
#ifndef SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_FIXEDWALL_H_
#define SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_FIXEDWALL_H_

#include "analysis/math/vector_distributed.h"
#include "esinfo/ecfinfo.h"
#include "feti/feti.h"

namespace espreso {

template <typename T>
struct FixedWall {

	void set(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet);
	void update(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet);

	std::vector<std::vector<esint> > cindex;
};

}




#endif /* SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_FIXEDWALL_H_ */
