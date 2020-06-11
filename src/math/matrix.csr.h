
#ifndef SRC_WRAPPERS_MATH_MATRIXCSR_H_
#define SRC_WRAPPERS_MATH_MATRIXCSR_H_

#include "matrix.h"
#include "data.matrix.csr.h"

namespace espreso {

class VectorDense;
class MatrixDense;
namespace MATH { struct CSRHandler; }

class MatrixCSR: public Matrix, public DataMatrixCSR
{
public:
	MatrixCSR();
	MatrixCSR(const MatrixCSR &other);
	MatrixCSR(const MatrixDense &other);
	MatrixCSR(esint nrows, esint ncols, esint nnz);
	MatrixCSR& operator=(const MatrixCSR &other);
	MatrixCSR& operator=(const MatrixDense &other);
	~MatrixCSR();

	MatrixCSR* copy();

	void structureUpdated();
	void swap(Matrix *other);
	void shallowCopy(const Matrix *other);
	void shallowCopyStructure(const Matrix *other);
	void deepCopy(const Matrix *other);
	void deepCopyStructure(const Matrix *other);
	void uniformCombination(const Matrix *first, const Matrix *second, int nfirst, int nsecond);

	void fill(double value);
	void fillData(const Matrix *in);
	void fillCombinedData(const Matrix *in,  esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum);

	void apply(const Vector *in, Vector *out);
	void apply(const Vectors *in, Vectors *out);

	void scale(double alpha);
	void add(double alpha, const Matrix *a);
	void sum(double alpha, const Matrix *a, double beta, const Matrix *b);
	void addToCombination(double scale, const Matrix *in,  esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum);
	void transpose();
	void transposeTo(MatrixCSR *out);

	void fillDiagonal(Vector *diagonal) const;

	double norm();

	void multiply(MatrixCSR &A, MatrixCSR &B, bool transposeA = false);
	void factorizeSymbolic();
	void factorizeNumeric();
	void factorize();
	void solve(const MatrixDense &rhs, MatrixDense &solution);
	void solve(const VectorDense &rhs, VectorDense &solution);
	void removeLower(MatrixType type);

	const char* name() const { return "MatrixCSR"; }

protected:
	MATH::CSRHandler *_math;
};

}

#endif /* SRC_WRAPPERS_MATH_MATRIXCSR_H_ */
