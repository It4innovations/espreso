
#ifndef SRC_WRAPPERS_MATH_VECTORSPARSEDISTRIBUTED_H_
#define SRC_WRAPPERS_MATH_VECTORSPARSEDISTRIBUTED_H_

#include "vector.h"
#include "data.vector.sparse.h"
#include "data.distributed.h"

namespace espreso {

class VectorSparseDistributed: public Vector, public DataVectorSparse, public DataDistributed
{
public:
	VectorSparseDistributed();
	VectorSparseDistributed(esint size, esint nnz, esint nhalo, esint nneighbors);

	~VectorSparseDistributed();

	VectorSparseDistributed* copy();

	void swap(Vector *other);
	void shallowCopy(const Vector *other);
	void shallowCopyStructure(const Vector *other);
	void shallowCopyFromHolder(const Vector *other, esint offset, esint nvectors);
	void deepCopy(const Vector *other);
	void deepCopyStructure(const Vector *other);
	void uniformCombination(const Vector *first, const Vector *second, int nfirst, int nsecond);

	void fill(double value);
	void fillData(const Vector *in);
	void fillCombinedValues(const Vector *in, esint offset, esint nsize, esint sumsize);
	void fillValuesFromCombination(const Vector *in, esint offset, esint nsize, esint sumsize);

	void scale(double alpha);
	void add(double alpha, const Vector *a);
	void sum(double alpha, const Vector *a, double beta, const Vector *b);
	void addToCombination(double alpha, const Vector *in, esint offset, esint nsize, esint sumsize);

	double norm();
	double max();
	double absmax();
	double dot(const Vector *other);

	void averageDuplications();

	const char* name() const { return "VectorSparseDistributed"; }
};

class VectorsSparseDistributed: public Vectors
{
public:
	VectorSparseDistributed* holder() { return reinterpret_cast<VectorSparseDistributed*>(_holder); }
	const VectorSparseDistributed* holder() const { return reinterpret_cast<VectorSparseDistributed*>(_holder); }

	VectorSparseDistributed& operator[](esint index) { return *reinterpret_cast<VectorSparseDistributed*>(_vectors[index]); }
	const VectorSparseDistributed& operator[](esint index) const { return *reinterpret_cast<VectorSparseDistributed*>(_vectors[index]); }

	VectorSparseDistributed* at(esint index) { return reinterpret_cast<VectorSparseDistributed*>(_vectors[index]); }
	const VectorSparseDistributed* at(esint index) const { return reinterpret_cast<VectorSparseDistributed*>(_vectors[index]); }

	VectorsSparseDistributed();
	~VectorsSparseDistributed();

	VectorsSparseDistributed* copy();

	void resize(esint nvectors, esint nnz, esint nhalo, esint nneighbors);
	void fillDistribution(esint nhalo, esint *halo, esint *distribution, int *neighbors);

	const char* name() const { return "VectorsSparseDistributed"; }

protected:
	Vector* create();
};

}



#endif /* SRC_WRAPPERS_MATH_VECTORSPARSEDISTRIBUTED_H_ */
