
#ifndef SRC_MATH2_UTILS_DECOMPOSED_DECOMPOSITION_H_
#define SRC_MATH2_UTILS_DECOMPOSED_DECOMPOSITION_H_

#include "basis/containers/serializededata.h"
#include "math2/utils/distributed/distribution.h"

namespace espreso {

struct DIndex { esint domain, index; };

inline bool operator==(const DIndex &left, const DIndex &right) { return left.domain == right.domain && left.index == right.index; }
inline bool operator!=(const DIndex &left, const DIndex &right) { return !(left == right); }
inline bool operator <(const DIndex &left, const DIndex &right) { return left.domain == right.domain ? left.index < right.index : left.domain < right.domain; }

struct DOFsDecomposition: public DOFsDistribution {
	serializededata<esint, DIndex> *dmap;
};

}

#endif /* SRC_MATH2_UTILS_DECOMPOSED_DECOMPOSITION_H_ */
