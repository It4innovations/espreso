#ifndef SRC_SOLVE_ACA_H_
#define SRC_SOLVE_ACA_H_

#include "math/math.h"
#include "esinfo/eslog.h"
#include "mem/morphing_system.h"

#include <vector>

#ifdef HAVE_MKL
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_lapack.h"
#include "mkl_spblas.h"
#endif

namespace espreso {
namespace MATH{
namespace SOLVER{

esint GMRESolverInternal_ACA(
	const MorphingMatrix &M,
	double *rhsVals, 
	double *results,
	double tolerance, 
	esint maxIterations, 
	esint &itercount
);
	
}//end of namespace SOLVER
}//end of namespace MATH
}//end of namespace espreso

#endif /* SRC_SOLVE_ACA_H_ */
