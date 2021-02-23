
#ifndef SRC_WRAPPERS_MATH_MATH_H_
#define SRC_WRAPPERS_MATH_MATH_H_

namespace espreso {

enum class MatrixType: int;

namespace MATH {

	class CSRHandlerData;
	class CSRHandler {
	public:
		CSRHandler(esint nrows, esint ncols, esint nnz, esint *rows, esint *cols, double *vals);
		~CSRHandler();

		void sizes(esint &rows, esint &cols, esint &nnz) const;
		void data(esint* &mRows, esint* &mCols, double* &mVals) const;

		CSRHandlerData *inner;
	};

	class IJVHandlerData;
	class IJVHandler {
	public:
		IJVHandler(esint nrows, esint ncols, esint nnz, esint *rows, esint *cols, double *vals);
		~IJVHandler();

		IJVHandlerData *inner;
	};

	void setNumberOfThreads(int numberOfThreads);

	void upCSRMatVecProduct(CSRHandler *A, double *vals, double *result);
	void upCSRMatVecProduct(esint rows, esint cols, esint *mRows, esint *mCols, float *mVals, float *vVals, float *result);
	void upCSRMatVecProduct(esint rows, esint cols, esint *mRows, esint *mCols, double *mVals, double *vVals, double *result);
	void CSRMatVecProduct(CSRHandler *A, double *vals, double *result);
	// C = op(A) * B with memory leak!!
	void CSRMatCSRMatProduct(CSRHandler *C, CSRHandler *A, CSRHandler *B, bool transposeA);

	// WARNING: do not use this function
	void CSRTranspose(CSRHandler *A, CSRHandler *At);

	void CSRTranspose(esint rows, esint cols, esint *aRows, esint *aCols, double *aVals, esint *bRows, esint *bCols, double *bVals);
	void CSRRemoveLower(esint rows, esint cols, esint *aRows, esint *aCols, double *aVals);

	void vecScale(esint size, float alpha, float *vVals);
	void vecScale(esint size, double alpha, double *vVals);

	double vecDot(esint size, double *vVals);
	double vecDot(esint size, double *a, double *b);
	void vecAdd(esint size, double *result, double alpha, double *other);
	void vecAddSparse(esint size, double *result, double alpha, esint *indices, double *other);
	void vecAddToSparse(esint size, double *result, double alpha, esint *indices, double *other);
	void vecSum(esint size, double *result, double alpha, double *a, double beta, double *b);

	void vecDenseToSparse(esint size, esint *indices, double *sparse, double *dense);

	double vecNorm(esint size, float *vVals);
	double vecNorm(esint size, double *vVals);

	void upDense3x3EigenValues(double *mVals, double *eigenValues);

	esint vecNormMaxIndex(esint size, float *vVals);
	esint vecNormMaxIndex(esint size, double *vVals);

	void DenseTranspose(esint rows, esint cols, double *vals);
	void DenseMinGeneralizedEigenVectors(esint msize, double *A, double *B, esint n, double *lambdas, double *vectors);

	// C = alpha * A * B + beta * C
	void DenseMatDenseMatRowMajorProduct(
			double alpha, bool transposeA, esint aRows, esint aCols, double* aVals,
			bool transposeB, esint bRows, esint bCols, double* bVals,
			double beta, double* cVals);
	// B = A \ B,  A in R^nrowsA.nrowsA, B in R^nrowsA.ncolsB
	void DenseMatDenseMatRowMajorSystemSolve(int nrowsA, int ncolsB, double *A, double *B);

	void CSRMatFactorizeSymbolic(MatrixType type, CSRHandler *A);
	void CSRMatFactorizeNumeric(MatrixType type, CSRHandler *A);
	void CSRMatFactorize(MatrixType type, CSRHandler *A);
	void CSRMatSolve(MatrixType type, CSRHandler *A, esint nrhs, double *rhs, double *solution);
	void CSRMatClearFactors(CSRHandler *A);

	inline double determinant2x2(double *values)
	{
		return values[0] * values[3] - values[1] * values[2];
	}

	inline double determinant3x3(double *values)
	{
		return
			+ values[0] * values[4] * values[8]
			+ values[1] * values[5] * values[6]
			+ values[2] * values[3] * values[7]
			- values[2] * values[4] * values[6]
			- values[1] * values[3] * values[8]
			- values[0] * values[5] * values[7];
	}

	inline void Dense2x2inverse(const double *m, double *inv, double det)
	{
		double detJx = 1 / det;
		inv[0] =   detJx * m[3];
		inv[1] = - detJx * m[1];
		inv[2] = - detJx * m[2];
		inv[3] =   detJx * m[0];
	}

	inline void Dense3x3inverse(const double *m, double *inv, double det)
	{
		double detJx = 1 / det;
		inv[0] = detJx * ( m[8] * m[4] - m[7] * m[5]);
		inv[1] = detJx * (-m[8] * m[1] + m[7] * m[2]);
		inv[2] = detJx * ( m[5] * m[1] - m[4] * m[2]);
		inv[3] = detJx * (-m[8] * m[3] + m[6] * m[5]);
		inv[4] = detJx * ( m[8] * m[0] - m[6] * m[2]);
		inv[5] = detJx * (-m[5] * m[0] + m[3] * m[2]);
		inv[6] = detJx * ( m[7] * m[3] - m[6] * m[4]);
		inv[7] = detJx * (-m[7] * m[0] + m[6] * m[1]);
		inv[8] = detJx * ( m[4] * m[0] - m[3] * m[1]);
	}

	namespace SOLVER {

		void GMRESUpCRSMat(
				esint rows, esint cols, esint *mRows, esint *mCols, double *mVals,
				double *rhsVals, double *results,
				double tolerance, esint maxIterations, esint &itercount);

		void GMRESDenseRowMajorMat(
				esint rows, esint cols, double *mVals,
				double *rhsVals, double *results,
				double tolerance, esint maxIterations, esint &itercount);

		void GMRESUpperSymetricColumnMajorMat(
				esint cols, double *mVals,
				double *rhsVals, double *results,
				double tolerance, esint maxIterations, esint &itercount);

		esint directUpperSymetricIndefiniteColumnMajor(
				esint cols, double *m_packed_values,
				esint nrhs, double *rhsVals);
	};
};

}



#endif /* SRC_WRAPPERS_MATH_MATH_H_ */
