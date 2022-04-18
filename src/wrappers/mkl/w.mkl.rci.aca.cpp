
#include "math/math.h"
#include "math/math.h"
#include "esinfo/eslog.h"
#include "morphing/morphing_system.h"

#include <vector>

#ifdef HAVE_MKL
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_lapack.h"
#include "mkl_spblas.h"
#endif

using namespace espreso;
//
//
//esint MATH::SOLVER::GMRESolverInternal_ACA(
//	const MorphingMatrix *M,
//	double *rhsVals,
//	double *results,
//	double tolerance,
//	esint maxIterations,
//	esint &itercount
//)
//{
//
//	esint niters = 0;
//	esint nToRestart = maxIterations;
//	double* mem_vec = new double[M->getNRows()];
//	// double* rhsVals_tmp = new double[M.getNRows()];
//
//
//	// M.applyPreconditioner(rhsVals, rhsVals_tmp, 1.0f, 0.0f, false);
//
//	// printf("c = [");
//	// for(int i = 0; i < M.getNRows(); ++i){
//	// 	printf("%f ", rhsVals_tmp[i]);
//	// }
//	// printf("\n];\n");
//
//#ifdef HAVE_MKL
//	//---------------------------------------------------------------------------
//	// Define arrays for the coefficient matrix
//	// Compressed sparse row storage is used for sparse representation
//	//---------------------------------------------------------------------------
//	esint size = 128;
//
//	esint cols = M->getNCols();
//
//	//---------------------------------------------------------------------------
//	// Allocate storage for the ?par parameters and the solution/rhs/residual vectors
//	//---------------------------------------------------------------------------
//	esint ipar[size];
//	double dpar[size];
//
//	ipar[14] = nToRestart;
//	std::vector<double> tmp(cols * (2 * ipar[14] + 1) + (ipar[14] * (ipar[14] + 9)) / 2 + 1);
//
//	//---------------------------------------------------------------------------
//	// Some additional variables to use with the RCI (P)FGMRES solver
//	//---------------------------------------------------------------------------
//	esint RCI_request, ivar;
//	ivar = cols;
//
//	double alpha = 1.0;
//	double beta = 0.0;
//
//	//---------------------------------------------------------------------------
//	// Initialize the initial guess
//	//---------------------------------------------------------------------------
//	for (esint i = 0; i < cols; i++) {
//		results[i] = 0.0;
//	}
//
//	//---------------------------------------------------------------------------
//	// Initialize the solver
//	//---------------------------------------------------------------------------
//	dfgmres_init(&ivar, results, rhsVals, &RCI_request, ipar, dpar, tmp.data());
//	// dfgmres_init(&ivar, results, rhsVals_tmp, &RCI_request, ipar, dpar, tmp.data());
//	if (RCI_request != 0) {
//		eslog::error("Something wrong happens while 'dfgmres_init' in ACA solver.\n");
//	}
//
//	//---------------------------------------------------------------------------
//	// Set the desired parameters:
//	// https://software.intel.com/en-us/node/521710
//	//---------------------------------------------------------------------------
//	ipar[4] = maxIterations;
//	dpar[0] = tolerance;
//
//	ipar[7] = 1;//1 -> performs test for maximal number of iterations: ipar[3] <= ipar[4]
//	ipar[8] = 1;//1 -> performs residual stopping test: dpar[4] <= dpar[3]
//	ipar[9] = 0;
//	ipar[10] = 0;//0 -> non-preconditioned version
//	ipar[11] = 1;//1 -> performs automatic test for zero norm of the currently generated vector: dpar[6] <= dpar[7], dpar[7] contains the tolerance value
//
//	//---------------------------------------------------------------------------
//	// Check the correctness and consistency of the newly set parameters
//	//---------------------------------------------------------------------------
//	dfgmres_check(&ivar, results, rhsVals, &RCI_request, ipar, dpar, tmp.data());
//	// dfgmres_check(&ivar, results, rhsVals_tmp, &RCI_request, ipar, dpar, tmp.data());
//	if (RCI_request != 0) {
//		eslog::error("Something wrong happens while 'dfgmres_check' in ACA solver.\n");
//	}
//
//	//---------------------------------------------------------------------------
//	// Compute the solution by RCI (P)FGMRES solver with preconditioning
//	// Reverse Communication starts here
//	//---------------------------------------------------------------------------
//	niters = 0;
//	while (true) {
//		dfgmres(&ivar, results, rhsVals, &RCI_request, ipar, dpar, tmp.data());
//		// dfgmres(&ivar, results, rhsVals_tmp, &RCI_request, ipar, dpar, tmp.data());
//
//		//---------------------------------------------------------------------------
//		// If RCI_request=0, then the solution was found with the required precision
//		//---------------------------------------------------------------------------
//		//std::cout<<"RCI "<<RCI_request<<std::endl;
//
//		if (RCI_request == 0) break;
//
//		//---------------------------------------------------------------------------
//		// If RCI_request=1, then compute the vector tmp[ipar[22]-1] = alpha*A*tmp[ipar[21]-1] + beta * tmp[ipar[22]-1]
//		//---------------------------------------------------------------------------
//		if (RCI_request == 1) {
//			M->apply(tmp.data() + ipar[21] - 1, tmp.data() + ipar[22] - 1, alpha, beta, false);
//			// M.apply(tmp.data() + ipar[21] - 1, mem_vec, 1.0f, 0.0f, false);
//			// M.applyPreconditioner(mem_vec, tmp.data() + ipar[22] - 1, alpha, beta, false);
//			niters++;
//			continue;
//		}
//
//
//		//---------------------------------------------------------------------------
//		// If RCI_request=3, then apply the preconditioner to tmp[ipar[21]-1] and store it into tmp[ipar[22]-1]
//		//---------------------------------------------------------------------------
//		if (RCI_request == 3) {
//			M->applyPreconditioner(tmp.data() + ipar[21] - 1, tmp.data() + ipar[22] - 1, 1.0f, 0.0f, false);
//			niters++;
//			continue;
//		}
//
//
//		break;
//	}
//	//---------------------------------------------------------------------------
//	// Reverse Communication ends here
//	// Get the current iteration number and the FGMRES solution (DO NOT FORGET to
//	// call dfgmres_get routine as computed_solution is still containing
//	// the initial guess!). Request to dfgmres_get to put the solution
//	// into vector computed_solution[N] via ipar[12]
//	//---------------------------------------------------------------------------
//	ipar[12] = 0;
//	dfgmres_get(&ivar, results, rhsVals, &RCI_request, ipar, dpar, tmp.data(), &itercount);
//	// dfgmres_get(&ivar, results, rhsVals_tmp, &RCI_request, ipar, dpar, tmp.data(), &itercount);
//	//---------------------------------------------------------------------------
//	// Print solution vector: computed_solution[N] and the number of iterations: itercount
//	//---------------------------------------------------------------------------
//	if (RCI_request != 0) {
//		eslog::error("Something wrong happens while 'dfgmres_get' in ACA solver.\n");
//	}
//#endif
//
//	delete [] mem_vec;
//	// delete [] rhsVals_tmp;
//
//	return niters;
//}
