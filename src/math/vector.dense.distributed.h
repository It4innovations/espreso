
#ifndef SRC_WRAPPERS_MATH_VECTORDENSEDISTRIBUTED_H_
#define SRC_WRAPPERS_MATH_VECTORDENSEDISTRIBUTED_H_

#include "vector.h"
#include "data.vector.dense.h"
#include "data.distributed.h"

namespace espreso {

class VectorFETI;

class VectorDenseDistributed: public Vector, public DataVectorDense, public DataDistributed
{
public:
	double& operator[](esint index) { return vals[index]; }
	const double& operator[](esint index) const { return vals[index]; }

	VectorDenseDistributed();
	VectorDenseDistributed(esint size, esint nhalo, esint nneighbors);

	~VectorDenseDistributed();

	VectorDenseDistributed* copy();

	void structeUpdated();
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

	void toFETI(VectorFETI *other) const;

	void gatherFromUpper();
	void scatterToUpper();

	const char* name() const { return "VectorDenseDistributed"; }
};

class VectorsDenseDistributed: public Vectors
{
public:
	VectorDenseDistributed* holder() { return reinterpret_cast<VectorDenseDistributed*>(_holder); }
	const VectorDenseDistributed* holder() const { return reinterpret_cast<VectorDenseDistributed*>(_holder); }

	VectorDenseDistributed& operator[](esint index) { return *reinterpret_cast<VectorDenseDistributed*>(_vectors[index]); }
	const VectorDenseDistributed& operator[](esint index) const { return *reinterpret_cast<VectorDenseDistributed*>(_vectors[index]); }

	VectorDenseDistributed* at(esint index) { return reinterpret_cast<VectorDenseDistributed*>(_vectors[index]); }
	const VectorDenseDistributed* at(esint index) const { return reinterpret_cast<VectorDenseDistributed*>(_vectors[index]); }

	VectorsDenseDistributed();
	~VectorsDenseDistributed();

	VectorsDenseDistributed* copy();

	void resize(esint size, esint nhalo, esint nneighbors);
	void structureUpdated();
	void fillDistribution(esint *halo, esint *distribution, int *neighbors);

	void gatherFromUpper();
	void scatterToUpper();

	const char* name() const { return "VectorsDenseDistributed"; }
protected:
	Vector* create();
};

}

#endif /* SRC_WRAPPERS_MATH_VECTORDENSEDISTRIBUTED_H_ */
