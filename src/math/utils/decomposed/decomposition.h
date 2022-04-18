
#ifndef SRC_MATH2_UTILS_DECOMPOSED_DECOMPOSITION_H_
#define SRC_MATH2_UTILS_DECOMPOSED_DECOMPOSITION_H_

#include "basis/containers/serializededata.h"
#include "math/utils/distributed/distribution.h"

namespace espreso {

struct DIndex { esint domain, index; };

inline bool operator==(const DIndex &left, const DIndex &right) { return left.domain == right.domain && left.index == right.index; }
inline bool operator!=(const DIndex &left, const DIndex &right) { return !(left == right); }
inline bool operator <(const DIndex &left, const DIndex &right) { return left.domain == right.domain ? left.index < right.index : left.domain < right.domain; }

struct DOFsDecomposition: public DOFsDistribution {
	esint dbegin, dend, dtotal;
	std::vector<esint> neighDomain; // first domain index per neighbor, the last is MY OFFSET

	std::vector<esint> sharedDOFs;
	serializededata<esint, DIndex> *dmap;

	bool ismy(esint domain) const { return dbegin <= domain && domain < dend; }
};

}

#endif /* SRC_MATH2_UTILS_DECOMPOSED_DECOMPOSITION_H_ */
