#include "solve_aca.h"

using namespace espreso;
// using namespace MATH;
// using namespace SOLVER;

esint MATH::SOLVER::GMRESolverInternal_ACA(
	const MorphingMatrix &M,
	double *rhsVals, 
	double *results,
	double tolerance, 
	esint maxIterations, 
	esint &itercount
)
{
	esint niters = 0;
#ifdef HAVE_MKL
	//---------------------------------------------------------------------------
	// Define arrays for the coefficient matrix
	// Compressed sparse row storage is used for sparse representation
	//---------------------------------------------------------------------------
	esint size = 128;
	
	esint rows = M.getNRows();
	esint cols = M.getNCols();

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

	double alpha = 1.0;
	double beta = 0.0;

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
		eslog::error("Something wrong happens while 'dfgmres_init' in ACA solver.\n");
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
		eslog::error("Something wrong happens while 'dfgmres_check' in ACA solver.\n");
	}

	//---------------------------------------------------------------------------
	// Compute the solution by RCI (P)FGMRES solver with preconditioning
	// Reverse Communication starts here
	//---------------------------------------------------------------------------
	niters = 0;
	while (true) {
		dfgmres(&ivar, results, rhsVals, &RCI_request, ipar, dpar, tmp.data());

		//---------------------------------------------------------------------------
		// If RCI_request=0, then the solution was found with the required precision
		//---------------------------------------------------------------------------
		//std::cout<<"RCI "<<RCI_request<<std::endl;

		if (RCI_request == 0) break;
		
		//---------------------------------------------------------------------------
		// If RCI_request=1, then compute the vector tmp[ipar[22]-1] = alpha*A*tmp[ipar[21]-1] + beta * tmp[ipar[22]-1]
		//---------------------------------------------------------------------------
		if (RCI_request == 1) {
			M.apply(tmp.data() + ipar[21] - 1, tmp.data() + ipar[22] - 1, alpha, beta, false);
			niters++;
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
		eslog::error("Something wrong happens while 'dfgmres_get' in ACA solver.\n");
	}
#endif

	return niters;
}
