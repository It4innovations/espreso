
#ifndef SRC_WRAPPERS_MATH_DATADECOMPOSITION_H_
#define SRC_WRAPPERS_MATH_DATADECOMPOSITION_H_

#include "domainindices.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

class DataDecompositionGather;
class VectorDenseFETI;

struct _DataDecomposition {
	// input
	esint rank, ranks, nneighbors;
	esint *distribution;
	int *neighbors;
	serializededata<esint, DI> *dmap;

	// computed
	esint doffset, nshared;
	esint *shared;

	// gathered data from all shared DOFs
	double *gathered;
	DataDecompositionGather *gather;

	_DataDecomposition();
	void alloc(esint ranks, esint nneighbors);
	void clear();
};

class DataDecomposition: public _DataDecomposition
{
public:
	enum class DUPLICATION {
		SPLIT,
		DUPLICATE,
		SPLIT_DOMAINS
	};

	DataDecomposition(DUPLICATION duplications);
	~DataDecomposition();

	void swap(DataDecomposition *other);
	void shallowCopy(const DataDecomposition *other);
	void shallowCopyStructure(const DataDecomposition *other);
	void deepCopy(const DataDecomposition *other);
	void deepCopyStructure(const DataDecomposition *other);
	void uniformCombination(const DataDecomposition *first, const DataDecomposition *second, int nfirst, int nsecond);
	void fillDecomposition(esint rank, esint ranks, esint nneighbors, esint *distribution, int *neighbors, const serializededata<esint, DI> *dmap);

	void allGather(const VectorDenseFETI &in) const;

	bool ismy(esint domain) const
	{
		return distribution[rank] <= domain && domain < distribution[rank + 1];
	}

	DUPLICATION duplications;
protected:
	serializededata<esint, DI>* combineDomainMap(const serializededata<esint, DI> *first, const serializededata<esint, DI> *second, int nfirst, int nsecond);
	void buildGatherData();

	_DataDecomposition _allocated;
};

}

#endif /* SRC_WRAPPERS_MATH_DATADECOMPOSITION_H_ */
