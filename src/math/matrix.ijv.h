
#ifndef SRC_WRAPPERS_MATH_MATRIXIJV_H_
#define SRC_WRAPPERS_MATH_MATRIXIJV_H_

#include "matrix.h"
#include "data.matrix.ijv.h"

namespace espreso {

class MatrixDense;
namespace MATH { struct IJVHandler; }

class MatrixIJV: public Matrix, public DataMatrixIJV
{
public:
	MatrixIJV();
	MatrixIJV(const MatrixIJV &other);
	MatrixIJV(const MatrixDense &other);
	MatrixIJV(esint nrows, esint ncols, esint nnz);
	MatrixIJV& operator=(const MatrixIJV &other);
	MatrixIJV& operator=(const MatrixDense &other);
	~MatrixIJV();

	MatrixIJV* copy();

	void structureUpdated();
	void swap(Matrix *other);
	void shallowCopy(const Matrix *other);
	void shallowCopyStructure(const Matrix *other);
	void deepCopy(const Matrix *other);
	void deepCopyStructure(const Matrix *other);
	void uniformCombination(const Matrix *first, const Matrix *second, int nfirst, int nsecond);

	void fill(double value);
	void fillData(const Matrix *in);
	void fillCombinedData(const Matrix *in, esint roffset, esint coffset, esint nsize, esint sumsize);

	void apply(const Vector *in, Vector *out);
	void apply(const Vectors *in, Vectors *out);

	void scale(double alpha);
	void add(double alpha, const Matrix *a);
	void sum(double alpha, const Matrix *a, double beta, const Matrix *b);
	void addToCombination(double scale, const Matrix *in, esint roffset, esint coffset, esint nsize, esint sumsize);
	void transpose();
	void transposeTo(MatrixIJV *out);

	void fillDiagonal(Vector *diagonal) const;

	double norm();

	void removeLower(MatrixType type);

	const char* name() const { return "MatrixIJV"; }

protected:
	MATH::IJVHandler *_math;
};

}

#endif /* SRC_WRAPPERS_MATH_MATRIXIJV_H_ */
