
#ifndef SOLVER_SPARSE_SPARSESOLVERS_H_
#define SOLVER_SPARSE_SPARSESOLVERS_H_

#if SOLVER == MKL
#include "cpu/mkl.h"
	typedef SparseSolverMKL SparseSolverCPU;
	typedef SparseSolverMKL SparseSolverAcc;


#elif SOLVER == PARDISO
#include "cpu/pardiso.h"
	typedef SparseSolverPardiso SparseSolverCPU;
	typedef SparseSolverPardiso SparseSolverAcc;


#elif SOLVER == MUMPS
#include "cpu/mumps.h"
	typedef SparseSolverMUMPS SparseSolverCPU;
	typedef SparseSolverMUMPS SparseSolverAcc;


#elif SOLVER == MIC
#include "acc/pardiso.h"
#include "acc/mic.h"
	typedef SparseSolverPARDISO SparseSolverCPU;
	typedef SparseSolverMIC SparseSolverAcc;


#elif SOLVER == CUDA
#include "acc/pardiso.h"
#include "acc/cuda.h"
	typedef SparseSolverPARDISO SparseSolverCPU;
	typedef SparseSolverCUDA SparseSolverAcc;


#else
#error "Incorrect user-supplied value for SOLVER."
#endif




#endif /* SOLVER_SPARSE_SPARSESOLVERS_H_ */
