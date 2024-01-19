
#ifndef SRC_ANALYSIS_LINEARSYSTEM_MATRICES_MATRIX_FETI_DECOMPOSITION_H_
#define SRC_ANALYSIS_LINEARSYSTEM_MATRICES_MATRIX_FETI_DECOMPOSITION_H_

#include "matrix_distributed.distribution.h"
#include "basis/containers/serializededata.h"

namespace espreso {

struct DIndex { esint domain, index; };

inline bool operator==(const DIndex &left, const DIndex &right) { return left.domain == right.domain && left.index == right.index; }
inline bool operator!=(const DIndex &left, const DIndex &right) { return !(left == right); }
inline bool operator <(const DIndex &left, const DIndex &right) { return left.domain == right.domain ? left.index < right.index : left.domain < right.domain; }

struct DOFsDecomposition: public DOFsDistribution {
	esint dbegin, dend, dtotal;
	std::vector<esint> neighDomain; // first domain index per neighbor, the last is MY OFFSET

	std::vector<esint> fixedDOFs;
	serializededata<esint, DIndex> *dmap = nullptr;

	bool ismy(esint domain) const
	{
		return dbegin <= domain && domain < dend;
	}

	int noffset(esint domain) const
	{
		size_t n = 0;
		while (n + 1 < neighbors.size() && neighDomain[n + 1] <= domain) { ++n; }
		return n;
	}

	~DOFsDecomposition()
	{
		if (dmap) delete dmap;
	}
};

}



#endif /* SRC_ANALYSIS_LINEARSYSTEM_MATRICES_MATRIX_FETI_DECOMPOSITION_H_ */
