
#ifndef SRC_WRAPPERS_MATH_MATRIX_H_
#define SRC_WRAPPERS_MATH_MATRIX_H_

#include "matrix.type.h"

namespace espreso {

class Vector;
class Vectors;

class Matrix
{
public:
	Matrix();
	Matrix(const Matrix &other);
	Matrix* shallowCopy();
	Matrix* shallowCopyStructure();
	Matrix* deepCopy();
	Matrix* deepCopyStructure();
	virtual ~Matrix();

	virtual Matrix* copy() =0;

	virtual void structureUpdated() =0;
	virtual void swap(Matrix *other) =0;
	virtual void shallowCopy(const Matrix *other) =0;
	virtual void shallowCopyStructure(const Matrix *other) =0;
	virtual void deepCopy(const Matrix *other) =0;
	virtual void deepCopyStructure(const Matrix *other) =0;
	virtual void uniformCombination(const Matrix *first, const Matrix *second, int nfirst, int nsecond) =0;

	virtual void fill(double value) =0;
	virtual void fillData(const Matrix *in) =0;
	virtual void fillCombinedData(const Matrix *in, esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum) =0;

	virtual void apply(const Vector *in, Vector *out) =0;
	virtual void apply(const Vectors *in, Vectors *out) =0;

	virtual void scale(double alpha) =0;
	virtual void add(double alpha, const Matrix *a) =0;
	virtual void sum(double alpha, const Matrix *a, double beta, const Matrix *b) =0;
	virtual void addToCombination(double scale, const Matrix *in, esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum) =0;

	virtual void fillDiagonal(Vector *diagonal) const =0;

	virtual double norm() =0;

	template <typename TMatrix>
	TMatrix* downcast()
	{
		if (dynamic_cast<TMatrix*>(this)) {
			return dynamic_cast<TMatrix*>(this);
		}
		downcastFailed(this, reinterpret_cast<Matrix*>(new TMatrix()));
		return 0;
	}

	template <typename TMatrix>
	const TMatrix* downcast() const
	{
		if (dynamic_cast<const TMatrix*>(this)) {
			return dynamic_cast<const TMatrix*>(this);
		}
		downcastFailed(this, reinterpret_cast<Matrix*>(new TMatrix()));
		return 0;
	}

	virtual const char* name() const =0;

	MatrixType type;
protected:
	void downcastFailed(const Matrix *m, const Matrix *target) const;

	void _assign(const Matrix *other);
	void _swap(Matrix *other);
};

}



#endif /* SRC_WRAPPERS_MATH_MATRIX_H_ */
