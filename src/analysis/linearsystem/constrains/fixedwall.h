
#ifndef SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_FIXEDWALL_H_
#define SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_FIXEDWALL_H_

#include "analysis/math/vector_distributed.h"
#include "esinfo/ecfinfo.h"
#include "feti/feti.h"

namespace espreso {

template <typename T>
struct FixedWall {

	static void set(const step::Step &step, FETI<T> &feti);
	static void update(const step::Step &step, FETI<T> &feti);
};

}




#endif /* SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_FIXEDWALL_H_ */
