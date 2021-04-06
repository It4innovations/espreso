
#ifndef SRC_WRAPPERS_MATH_MATRIXDENSE_H_
#define SRC_WRAPPERS_MATH_MATRIXDENSE_H_

#include "matrix.h"
#include "data.matrix.dense.h"

namespace espreso {

class MatrixCSR;

class MatrixDense: public Matrix, public DataMatrixDense
{
public:
	class MatrixDenseRow
	{
	public:
		MatrixDenseRow(double *vals): vals(vals) {}

		double& operator[](esint col) { return vals[col]; }
		const double& operator[](esint col) const { return vals[col]; }
	protected:
		double *vals;
	};

	MatrixDenseRow operator[](esint row) { return MatrixDenseRow(vals + ncols * row); }
	const MatrixDenseRow operator[](esint row) const { return MatrixDenseRow(vals + ncols * row); }

	double& operator()(esint row, esint col) { return vals[ncols * row + col]; }
	const double& operator()(esint row, esint col) const { return vals[ncols * row + col]; }

	double& at(esint row, esint col) { return vals[ncols * row + col]; }
	const double& at(esint row, esint col) const { return vals[ncols * row + col]; }

	MatrixDense();
	MatrixDense(esint nrows, esint ncols);
	MatrixDense(const MatrixDense &other);
	MatrixDense(const MatrixCSR &other);
	MatrixDense& operator=(const MatrixDense &other);
	MatrixDense& operator=(const MatrixCSR &other);
	~MatrixDense();

	MatrixDense* copy();

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
	void addToCombination(double alpha, const Matrix *in, esint roffset, esint coffset, esint nsize, esint sumsize);

	void fillDiagonal(Vector *diagonal) const;

	double norm();
	esint nnz(double eps) const;

	void multiply(
			const MatrixDense &A, const MatrixDense &B,
			double alpha = 1, double beta = 0,
			bool transposeA = false, bool transposeB = false);

	void multiply(
			const MatrixDense &A,
			esint brows, esint bcols, double* bvals,
			double alpha = 1, double beta = 0,
			bool transposeA = false, bool transposeB = false);

	void multiply(
			double *A,
			esint arows, esint acols, esint brows, esint bcols, double* bvals,
			double alpha = 1, double beta = 0,
			bool transposeA = false, bool transposeB = false);

	void systemSolve(const MatrixDense &A, const MatrixDense &B);

	const char* name() const { return "MatrixDense"; }
};

}



#endif /* SRC_WRAPPERS_MATH_MATRIXDENSE_H_ */
