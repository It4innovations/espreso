
#ifndef SRC_MATH2_FETI_LMAP_H_
#define SRC_MATH2_FETI_LMAP_H_

#include <vector>

namespace espreso {

struct LMAP {
	enum: int {
		LOCAL = -1,
		DIRICHLET = -2
	};

	int from, to, offset, neigh;
};


}

#endif /* SRC_MATH2_FETI_LMAP_H_ */
