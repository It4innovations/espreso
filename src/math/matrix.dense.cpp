
#include "matrix.dense.h"
#include "matrix.csr.h"
#include "vector.dense.h"
#include "math.h"
#include "esinfo/eslog.h"

using namespace espreso;

MatrixDense::MatrixDense()
{

}

MatrixDense::MatrixDense(esint nrows, esint ncols)
: DataMatrixDense(nrows, ncols)
{

}

MatrixDense::MatrixDense(esint nrows, esint ncols, double *vals)
: DataMatrixDense(nrows, ncols, vals)
{

}

MatrixDense::MatrixDense(const MatrixDense &other)
: Matrix(other), DataMatrixDense(other)
{

}

MatrixDense::MatrixDense(const MatrixCSR &other)
: Matrix(other)
{
	resize(other.nrows, other.ncols);
	fill(0);
	esint indexing = other.rows[0];
	for (esint r = 0; r < other.nrows; ++r) {
		for (esint c = other.rows[r]; c < other.rows[r + 1]; ++c) {
			at(r, other.cols[c - indexing] - indexing) = other.vals[c - indexing];
		}
	}
}

MatrixDense& MatrixDense::operator=(const MatrixDense &other)
{
	if (this != &other) {
		deepCopy(&other);
	}
	return *this;
}

MatrixDense& MatrixDense::operator=(const MatrixCSR &other)
{
	Matrix::_assign(&other);
	resize(other.nrows, other.ncols);
	fill(0);
	esint indexing = other.rows[0];
	for (esint r = 0; r < other.nrows; ++r) {
		for (esint c = other.rows[r]; c < other.rows[r + 1]; ++c) {
			at(r, other.cols[c - indexing] - indexing) = other.vals[c - indexing];
		}
	}
	return *this;
}

MatrixDense::~MatrixDense()
{

}

MatrixDense* MatrixDense::copy()
{
	return new MatrixDense();
}

void MatrixDense::set(esint nrows, esint ncols, double *vals)
{
	DataMatrixDense::set(nrows, ncols, vals);
}

void MatrixDense::structureUpdated()
{

}

void MatrixDense::swap(Matrix *other)
{
	MatrixDense *_other = other->downcast<MatrixDense>();
	Matrix::_swap(other);
	DataMatrixDense::swap(_other);
}

void MatrixDense::shallowCopy(const Matrix *other)
{
	Matrix::_assign(other);
	DataMatrixDense::shallowCopy(other->downcast<MatrixDense>());
}

void MatrixDense::shallowCopyStructure(const Matrix *other)
{
	deepCopy(other); // there is not any structure
}

void MatrixDense::deepCopy(const Matrix *other)
{
	Matrix::_assign(other);
	DataMatrixDense::deepCopy(other->downcast<MatrixDense>());
}

void MatrixDense::deepCopyStructure(const Matrix *other)
{
	deepCopy(other);
}

void MatrixDense::uniformCombination(const Matrix *first, const Matrix *second, int nfirst, int nsecond)
{
	DataMatrixDense::uniformCombination(first->downcast<MatrixDense>(), second->downcast<MatrixDense>(), nfirst, nsecond);
}

void MatrixDense::fill(double value)
{
	DataMatrixDense::fill(value);
}

void MatrixDense::fillData(const Matrix *in)
{
	const MatrixDense *_in = in->downcast<MatrixDense>();
	DataMatrixDense::fillValues(_in->vals);
}

void MatrixDense::fillCombinedData(const Matrix *in,  esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum)
{
	const MatrixDense *_in = in->downcast<MatrixDense>();
	DataMatrixDense::fillCombinedValues(_in, roffset, coffset, rsize, csize, rsum, csum);
}

void MatrixDense::apply(const Vector *in, Vector *out)
{
	eslog::internalFailure("call empty function.\n");
}

void MatrixDense::apply(const Vectors *in, Vectors *out)
{
	eslog::internalFailure("call empty function.\n");
}

void MatrixDense::scale(double alpha)
{
	MATH::vecScale(nrows * ncols, alpha, vals);
}

void MatrixDense::add(double alpha, const Matrix *a)
{
	const MatrixDense *_a = a->downcast<MatrixDense>();
	MATH::vecAdd(nrows * ncols, vals, alpha, _a->vals);
}

void MatrixDense::sum(double alpha, const Matrix *a, double beta, const Matrix *b)
{
	const MatrixDense *_a = a->downcast<MatrixDense>();
	const MatrixDense *_b = b->downcast<MatrixDense>();
	MATH::vecSum(nrows * ncols, vals, alpha, _a->vals, beta, _b->vals);
}

void MatrixDense::addToCombination(double alpha, const Matrix *in,  esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum)
{
	DataMatrixDense::addToCombination(alpha, in->downcast<MatrixDense>(), roffset, coffset, rsize, csize, rsum, csum);
}

void MatrixDense::transpose()
{
	MATH::DenseTranspose(nrows, ncols, vals);
}

void MatrixDense::minGeneralizedEigenValues(double *B, esint n, double *lambdas, double *vectors)
{
	if (nrows == ncols) {
		MATH::DenseMinGeneralizedEigenVectors(nrows, vals, B, n, lambdas, vectors);
	} else {
		eslog::failure("Cannot compute generalized eigen problem of non-square matrices\n");
	}
}

void MatrixDense::fillDiagonal(Vector *diagonal) const
{
	MatrixDense::fillDiagonal(diagonal->downcast<VectorDense>());
}

double MatrixDense::norm()
{
	return MATH::vecNorm(nrows * ncols, vals);
}

esint MatrixDense::nnz(double eps) const
{
	esint nnz = 0;
	for (esint r = 0; r < nrows; r++) {
		for (esint c = 0; c < ncols; c++) {
			if (at(r, c) < -eps || eps < at(r, c)) {
				++nnz;
			}
		}
	}
	return nnz;
}

void MatrixDense::multiply(
		const MatrixDense &A, const MatrixDense &B,
		double alpha, double beta,
		bool transposeA, bool transposeB)
{
	resize(transposeA ? A.ncols : A.nrows, transposeB ? B.nrows : B.ncols);
	MATH::DenseMatDenseMatRowMajorProduct(
			alpha, transposeA, A.nrows, A.ncols, A.vals,
			transposeB, B.nrows, B.ncols, B.vals,
			beta, vals);
}

void MatrixDense::multiply(
		const MatrixDense &A,
		esint brows, esint bcols, double* bvals,
		double alpha, double beta,
		bool transposeA, bool transposeB)
{
	resize(transposeA ? A.ncols : A.nrows, transposeB ? brows : bcols);
	MATH::DenseMatDenseMatRowMajorProduct(
			alpha, transposeA, A.nrows, A.ncols, A.vals,
			transposeB, brows, bcols, bvals,
			beta, vals);
}

void MatrixDense::multiply(
		double *A,
		esint arows, esint acols, esint brows, esint bcols, double* bvals,
		double alpha, double beta,
		bool transposeA, bool transposeB)
{
	resize(transposeA ? acols : arows, transposeB ? brows : bcols);
	MATH::DenseMatDenseMatRowMajorProduct(
			alpha, transposeA, arows, acols, A,
			transposeB, brows, bcols, bvals,
			beta, vals);
}

void MatrixDense::systemSolve(const MatrixDense &A, const MatrixDense &B) {
	MatrixDense AA = A;
	deepCopy(&B);
	vals = B.vals;
	MATH::DenseMatDenseMatRowMajorSystemSolve(A.nrows, B.ncols, AA.vals, vals);
}



