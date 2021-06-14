
#ifndef SRC_WRAPPERS_MATH_VECTORFETI_H_
#define SRC_WRAPPERS_MATH_VECTORFETI_H_

#include "data.decomposition.h"
#include "vector.dense.h"

namespace espreso {

class VectorFETI: public Vector, public DataDecomposition
{
public:
	VectorFETI();
	virtual ~VectorFETI();

	Vector& operator[](esint domain) { return *vectors[domain]; }
	const Vector& operator[](esint domain) const { return *vectors[domain]; }

	Vector* at(esint domain) { return vectors[domain]; }
	const Vector* at(esint domain) const { return vectors[domain]; }

	virtual VectorFETI* copy();

	virtual void setDuplications(DataDecomposition::DUPLICATION duplications);
	virtual void initDomains(DataDecomposition::DUPLICATION duplications, esint domains);

	virtual void swap(Vector *other);
	virtual void shallowCopy(const Vector *other);
	virtual void shallowCopyStructure(const Vector *other);
	virtual void shallowCopyFromHolder(const Vector *other, esint offset, esint nvectors);
	virtual void deepCopy(const Vector *other);
	virtual void deepCopyStructure(const Vector *other);
	virtual void uniformCombination(const Vector *first, const Vector *second, int nfirst, int nsecond);

	virtual void fill(double value);
	virtual void fillData(const Vector *in);
	virtual void fillCombinedValues(const Vector *in, esint offset, esint nsize, esint sumsize);
	virtual void fillValuesFromCombination(const Vector *in, esint offset, esint nsize, esint sumsize);

	virtual void scale(double alpha);
	virtual void add(double alpha, const Vector *a);
	virtual void sum(double alpha, const Vector *a, double beta, const Vector *b);
	virtual void addToCombination(double alpha, const Vector *in, esint offset, esint nsize, esint sumsize);

	virtual double norm();
	virtual double max();
	virtual double absmax();
	virtual double dot(const Vector *other);

	virtual void fromFETI(VectorFETI *other) const;

	esint domains;
	Vector **vectors;

	const char* name() const { return "VectorFETI"; }

protected:
	virtual Vector* create();
};

class VectorsFETI: public Vectors
{
public:
	// used by print
	struct Domain {
		Domain(const VectorsFETI* data, esint domain): data(data), domain(domain) {}
		const VectorsFETI* data;
		esint domain;
	};

	Vector& operator[](esint index) { return *_vectors[index]; }
	const Vector& operator[](esint index) const { return *_vectors[index]; }

	Vector* at(esint index) { return _vectors[index]; }
	const Vector* at(esint index) const { return _vectors[index]; }

	VectorsFETI();
	virtual ~VectorsFETI();

	virtual VectorsFETI* copy();

	virtual void setDuplications(DataDecomposition::DUPLICATION duplications);

	const char* name() const { return "VectorsFETI"; }
protected:
	virtual Vector* create();
};

}

#endif /* SRC_WRAPPERS_MATH_VECTORFETI_H_ */
