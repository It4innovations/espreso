
#ifndef SRC_WRAPPERS_MATH_VECTOR_H_
#define SRC_WRAPPERS_MATH_VECTOR_H_

namespace espreso {

class Vector
{
public:
	virtual ~Vector();
	Vector* shallowCopy();
	Vector* shallowCopyStructure();
	Vector* deepCopy();
	Vector* deepCopyStructure();

	virtual Vector* copy() =0;

	virtual void swap(Vector *other) =0;
	virtual void shallowCopy(const Vector *other) =0;
	virtual void shallowCopyStructure(const Vector *other) =0;
	virtual void shallowCopyFromHolder(const Vector *other, esint offset, esint nvectors) =0;
	virtual void deepCopy(const Vector *other) =0;
	virtual void deepCopyStructure(const Vector *other) =0;
	virtual void uniformCombination(const Vector *first, const Vector *second, int nfirst, int nsecond) =0;

	virtual void fill(double value) =0;
	virtual void fillData(const Vector *in) =0;
	virtual void fillCombinedValues(const Vector *in, esint offset, esint nsize, esint sumsize) =0;
	virtual void fillValuesFromCombination(const Vector *in, esint offset, esint nsize, esint sumsize) =0;

	virtual void scale(double alpha) =0;
	virtual void add(double alpha, const Vector *a) =0;
	virtual void sum(double alpha, const Vector *a, double beta, const Vector *b) =0;
	virtual void addToCombination(double alpha, const Vector *in, esint offset, esint nsize, esint sumsize) =0;

	virtual double norm() =0;
	virtual double max() =0;
	virtual double absmax() =0;
	virtual double dot(const Vector *other) =0;

	template <typename TVector>
	TVector* downcast()
	{
		if (dynamic_cast<TVector*>(this)) {
			return dynamic_cast<TVector*>(this);
		}
		downcastFailed(this, reinterpret_cast<Vector*>(new TVector()));
		return 0;
	}

	template <typename TVector>
	const TVector* downcast() const
	{
		if (dynamic_cast<const TVector*>(this)) {
			return dynamic_cast<const TVector*>(this);
		}
		downcastFailed(this, reinterpret_cast<Vector*>(new TVector()));
		return 0;
	}

	virtual const char* name() const =0;
protected:
	void downcastFailed(const Vector *v, const Vector *target) const;
};

class Vectors
{
public:
	Vector& operator[](esint index) { return *_vectors[index]; }
	const Vector& operator[](esint index) const { return *_vectors[index]; }

	Vector* at(esint index) { return _vectors[index]; }
	const Vector* at(esint index) const { return _vectors[index]; }

	Vectors();
	virtual ~Vectors();

	Vectors* shallowCopy();
	Vectors* shallowCopyStructure();
	Vectors* deepCopy();
	Vectors* deepCopyStructure();

	virtual Vectors* copy() =0;

	virtual void initVectors(esint nvectors);

	virtual void swap(Vectors *other);
	virtual void shallowCopy(const Vectors *other);
	virtual void shallowCopyStructure(const Vectors *other);
	virtual void deepCopy(const Vectors *other);
	virtual void deepCopyStructure(const Vectors *other);
	virtual void uniformCombination(const Vector *first, const Vector *second, int nfirst, int nsecond);
	virtual void uniformCombination(const Vectors *first, const Vectors *second, int nfirst, int nsecond);

	virtual void fill(double value);
	virtual void fillData(const Vectors *in);
	virtual void fillCombinedValues(const Vectors *in, esint offset, esint nsize, esint sumsize);
	virtual void fillValuesFromCombination(const Vectors *in, esint offset, esint nsize, esint sumsize);

	virtual void scale(double alpha);
	virtual void add(double alpha, const Vectors *a);
	virtual void sum(double alpha, const Vectors *a, double beta, const Vectors *b);
	virtual void addToCombination(double alpha, const Vectors *in, esint offset, esint nsize, esint sumsize);

	template <typename TVectors>
	TVectors* downcast()
	{
		if (dynamic_cast<TVectors*>(this)) {
			return dynamic_cast<TVectors*>(this);
		}
		downcastFailed(this, reinterpret_cast<Vectors*>(new TVectors()));
		return 0;
	}

	template <typename TVectors>
	const TVectors* downcast() const
	{
		if (dynamic_cast<const TVectors*>(this)) {
			return dynamic_cast<const TVectors*>(this);
		}
		downcastFailed(this, reinterpret_cast<Vectors*>(new TVectors()));
		return 0;
	}

	virtual const char* name() const =0;

	esint nvectors;
protected:
	Vector *_holder;
	Vector **_vectors;

	virtual Vector* create() =0;
	void downcastFailed(const Vectors *v, const Vectors *target) const;
};

}



#endif /* SRC_WRAPPERS_MATH_VECTOR_H_ */
