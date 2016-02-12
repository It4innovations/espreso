
#ifndef SOLVER_SPECIFIC_SPARSESOLVERS_H_
#define SOLVER_SPECIFIC_SPARSESOLVERS_H_


#if defined(SOLVER_MKL)
#include "cpu/solvermkl.h"
	typedef SparseSolverMKL SparseSolverCPU;
	typedef SparseSolverMKL SparseSolverAcc;


#elif defined(SOLVER_PARDISO)
#include "cpu/solverpardiso.h"
	typedef SparseSolverPardiso SparseSolverCPU;
	typedef SparseSolverPardiso SparseSolverAcc;

#elif defined(SOLVER_MUMPS)
#include "cpu/solvermumps.h"
	typedef SparseSolverMUMPS SparseSolverCPU;
	typedef SparseSolverMUMPS SparseSolverAcc;

#elif defined(SOLVER_MIC)
#include "cpu/solverpardiso.h"
#include "acc/mic.h"
	typedef SparseSolverPARDISO SparseSolverCPU;
	typedef SparseSolverMIC SparseSolverAcc;


#elif defined(SOLVER_CUDA)
#include "cpu/solverpardiso.h"
#include "acc/cuda.h"
	typedef SparseSolverPARDISO SparseSolverCPU;
	typedef SparseSolverCUDA SparseSolverAcc;


#else
#error "Incorrect user-supplied value for SOLVER. Check your build.config script."
#endif




#endif /* SOLVER_SPECIFIC_SPARSESOLVERS_H_ */
