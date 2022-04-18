
#ifndef SRC_WRAPPERS_MATH_VECTORSPARSE_H_
#define SRC_WRAPPERS_MATH_VECTORSPARSE_H_

#include "vector.h"
#include "data.vector.sparse.h"

namespace espreso {

class VectorDense;
class VectorsDense;
class VectorFETI;

class VectorSparse: public Vector, public DataVectorSparse
{
public:
	VectorSparse();
	VectorSparse(esint size, esint nnz);

	~VectorSparse();

	VectorSparse* copy();

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
	void toCombinedFETI(VectorFETI *other, esint offset, esint nsize, esint sumsize) const;

	const char* name() const { return "VectorSparse"; }
};

class VectorsSparse: public Vectors
{
public:
	VectorSparse* holder() { return reinterpret_cast<VectorSparse*>(_holder); }
	const VectorSparse* holder() const { return reinterpret_cast<VectorSparse*>(_holder); }

	VectorSparse& operator[](esint index) { return *reinterpret_cast<VectorSparse*>(_vectors[index]); }
	const VectorSparse& operator[](esint index) const { return *reinterpret_cast<VectorSparse*>(_vectors[index]); }

	VectorSparse* at(esint index) { return reinterpret_cast<VectorSparse*>(_vectors[index]); }
	const VectorSparse* at(esint index) const { return reinterpret_cast<VectorSparse*>(_vectors[index]); }

	VectorsSparse();
	~VectorsSparse();

	VectorsSparse* copy();

	void resize(esint size, esint nnz);

	const char* name() const { return "VectorsSparse"; }

protected:
	Vector* create();
};

}


#endif /* SRC_WRAPPERS_MATH_VECTORSPARSE_H_ */
