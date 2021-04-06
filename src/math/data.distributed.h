
#ifndef SRC_WRAPPERS_MATH_DATADISTRIBUTED_H_
#define SRC_WRAPPERS_MATH_DATADISTRIBUTED_H_

namespace espreso {

struct DataSynchronization;
struct DataMV;

struct _DataDistributed {
	esint nhalo, ranks, nneighbors;
	esint *halo, *distribution, *nintervals;
	int *neighbors;

	DataSynchronization *sync;
	DataMV *mv;

	_DataDistributed();
	void alloc(esint nhalo, esint ranks, esint nneighbors);
	void clear();
};

class DataDistributed: public _DataDistributed
{
public:
	DataDistributed();
	DataDistributed(esint nhalo, esint ranks, esint nneighbors);
	~DataDistributed();

	void resize(esint nhalo, esint ranks, esint nneighbors);
	void swap(DataDistributed *other);
	void shallowCopy(const DataDistributed *other);
	void deepCopy(const DataDistributed *other);
	void uniformCombination(const DataDistributed *first, const DataDistributed *second, esint nfirst, esint nsecond);

	void fillDistribution(esint *halo, esint *distribution, int *neighbors);
protected:
	_DataDistributed _allocated;
};

}

#endif /* SRC_WRAPPERS_MATH_DATADISTRIBUTED_H_ */
