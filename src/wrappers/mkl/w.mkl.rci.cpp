
#include "math/math.h"
#include "math2/math2.h"
#include "esinfo/eslog.h"

#include <vector>

#ifdef HAVE_MKL
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_lapack.h"
#include "mkl_spblas.h"
#endif

using namespace espreso;

enum class SOLVER_INTERNAL_TYPE {
	DCSRGEMV,
	DGEMV,
	DSPMV,
};

void GMRESolverInternal(SOLVER_INTERNAL_TYPE type,
						esint rows, esint cols, esint *mRows, esint *mCols, double *mVals,
						double *rhsVals, double *results,
						double tolerance, esint maxIterations, esint &itercount)
{
#ifdef HAVE_MKL
	sparse_matrix_t inspector;
	matrix_descr descr;
	if (type == SOLVER_INTERNAL_TYPE::DCSRGEMV) {
		mkl_sparse_d_create_csr(&inspector, SPARSE_INDEX_BASE_ONE, rows, cols, mRows, mRows + 1, mCols, mVals);
		descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
		descr.mode = SPARSE_FILL_MODE_UPPER;
	}

	//---------------------------------------------------------------------------
	// Define arrays for the coefficient matrix
	// Compressed sparse row storage is used for sparse representation
	//---------------------------------------------------------------------------
	esint size = 128;

	//---------------------------------------------------------------------------
	// Allocate storage for the ?par parameters and the solution/rhs/residual vectors
	//---------------------------------------------------------------------------
	esint ipar[size];
	double dpar[size];
	std::vector<double> tmp(cols * (2 * cols + 1) + (cols * (cols + 9)) / 2 + 1);

	//---------------------------------------------------------------------------
	// Some additional variables to use with the RCI (P)FGMRES solver
	//---------------------------------------------------------------------------
	esint RCI_request, ivar;
	ivar = cols;

	char uplo = 'U';
	char transa = 'N';

	double alpha = 1.0;
	double beta = 0.0;

	esint incx = 1;
	esint incy = 1;

	//---------------------------------------------------------------------------
	// Initialize the initial guess
	//---------------------------------------------------------------------------
	for (esint i = 0; i < cols; i++) {
		results[i] = 0.0;
	}

	//---------------------------------------------------------------------------
	// Initialize the solver
	//---------------------------------------------------------------------------
	dfgmres_init(&ivar, results, rhsVals, &RCI_request, ipar, dpar, tmp.data());
	if (RCI_request != 0) {
		eslog::error("Something wrong happens while 'dfgmres_init' in solver.\n");
	}

	//---------------------------------------------------------------------------
	// Set the desired parameters:
	// https://software.intel.com/en-us/node/521710
	//---------------------------------------------------------------------------
	ipar[4] = maxIterations;
	dpar[0] = tolerance;

	ipar[7] = 1;
	ipar[8] = 1;
	ipar[9] = 0;
	ipar[10] = 0;
	ipar[11] = 1;

	//---------------------------------------------------------------------------
	// Check the correctness and consistency of the newly set parameters
	//---------------------------------------------------------------------------
	dfgmres_check(&ivar, results, rhsVals, &RCI_request, ipar, dpar, tmp.data());
	if (RCI_request != 0) {
		eslog::error("Something wrong happens while 'dfgmres_check' in solver.\n");
	}

	//---------------------------------------------------------------------------
	// Compute the solution by RCI (P)FGMRES solver with preconditioning
	// Reverse Communication starts here
	//---------------------------------------------------------------------------
	while (true) {
		dfgmres(&ivar, results, rhsVals, &RCI_request, ipar, dpar, tmp.data());
		//---------------------------------------------------------------------------
		// If RCI_request=0, then the solution was found with the required precision
		//---------------------------------------------------------------------------
		//std::cout<<"RCI "<<RCI_request<<std::endl;

		if (RCI_request == 0) break;
		//---------------------------------------------------------------------------
		// If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
		// and put the result in vector tmp[ipar[22]-1]
		//---------------------------------------------------------------------------
		if (RCI_request == 1) {
			switch (type) {
				case SOLVER_INTERNAL_TYPE::DSPMV: {
					dspmv(&uplo, &cols, &alpha, mVals, tmp.data() + ipar[21] - 1, &incx, &beta, tmp.data() + ipar[22] - 1, &incy);
				} break;

				case SOLVER_INTERNAL_TYPE::DGEMV: {
					dgemv(&transa, &rows, &cols, &alpha, mVals, &rows, tmp.data() + ipar[21] - 1, &incx, &beta, tmp.data() + ipar[22] - 1, &incy );
				}

				case SOLVER_INTERNAL_TYPE::DCSRGEMV: {
					if (transa == 'T') {
						mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, 1, inspector, descr, tmp.data() + ipar[21] - 1, 0, tmp.data() + ipar[22] - 1);
					} else {
						mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, inspector, descr, tmp.data() + ipar[21] - 1, 0, tmp.data() + ipar[22] - 1);
					}
				}

				default: eslog::internalFailure("no such internal solver.\n");
			}
			continue;
		}

		break;
	}
	//---------------------------------------------------------------------------
	// Reverse Communication ends here
	// Get the current iteration number and the FGMRES solution (DO NOT FORGET to
	// call dfgmres_get routine as computed_solution is still containing
	// the initial guess!). Request to dfgmres_get to put the solution
	// into vector computed_solution[N] via ipar[12]
	//---------------------------------------------------------------------------
	ipar[12] = 0;
	dfgmres_get(&ivar, results, rhsVals, &RCI_request, ipar, dpar, tmp.data(), &itercount);
	//---------------------------------------------------------------------------
	// Print solution vector: computed_solution[N] and the number of iterations: itercount
	//---------------------------------------------------------------------------
	if (RCI_request != 0) {
		eslog::error("Something wrong happens while 'dfgmres_get' in solver.\n");
	}
	if (type == SOLVER_INTERNAL_TYPE::DCSRGEMV) {
		mkl_sparse_destroy(inspector);
	}
#endif
}

void MATH::SOLVER::GMRESUpCRSMat(
		esint rows, esint cols, esint *mRows, esint *mCols, double *mVals,
		double *rhsVals, double *results,
		double tolerance, esint maxIterations, esint &itercount)
{
#ifdef HAVE_MKL
	GMRESolverInternal(SOLVER_INTERNAL_TYPE::DCSRGEMV,
			rows, cols, mRows, mCols, mVals,
			rhsVals, results,
			tolerance, maxIterations, itercount);
#endif
}

void MATH::SOLVER::GMRESDenseRowMajorMat(
		esint rows, esint cols, double *mVals,
		double *rhsVals, double *results,
		double tolerance, esint maxIterations, esint &itercount)
{
#ifdef HAVE_MKL
	GMRESolverInternal(SOLVER_INTERNAL_TYPE::DGEMV,
							rows, cols, NULL, NULL, mVals,
							rhsVals, results,
							tolerance, maxIterations, itercount);
#endif
}


void MATH::SOLVER::GMRESUpperSymetricColumnMajorMat(
		esint cols, double *mVals,
		double *rhsVals, double *results,
		double tolerance, esint maxIterations, esint &itercount)
{
#ifdef HAVE_MKL
	GMRESolverInternal(SOLVER_INTERNAL_TYPE::DSPMV,
							cols, cols, NULL, NULL, mVals,
							rhsVals, results,
							tolerance, maxIterations, itercount);
#endif
}

esint MATH::SOLVER::directUpperSymetricIndefiniteColumnMajor(
		esint cols, double *m_packed_values,
		esint nrhs, double *rhsVals)
{
	esint info = 0;
#ifdef HAVE_MKL
	char U = 'U';
	std::vector<esint>  m_ipiv(cols);
	dsptrf(&U, &cols, &m_packed_values[0], &m_ipiv[0], &info);
	if (info == 0) {
		dsptrs(&U, &cols, &nrhs, &m_packed_values[0], &m_ipiv[0], &rhsVals[0], &cols, &info);
	}
#endif
	return info;
}


