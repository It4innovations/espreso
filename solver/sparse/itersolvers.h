
#ifndef SOLVER_SPARSE_ITERSOLVERS_H_
#define SOLVER_SPARSE_ITERSOLVERS_H_

#if defined(SOLVER_MKL)
#include "cpu/itersolvercpu.h"
	typedef IterSolverCPU	IterSolver;


#elif defined(SOLVER_PARDISO)
#include "cpu/itersolvercpu.h"
	typedef IterSolverCPU	IterSolver;

#elif defined(SOLVER_MUMPS)
#include "cpu/itersolvercpu.h"
	typedef IterSolverCPU	IterSolver;

#elif defined(SOLVER_MIC)
#include "acc/itersolveracc.h"
	typedef IterSolverAcc	IterSolver;

#elif defined(SOLVER_CUDA)
#include "acc/itersolveracc.h"
	typedef IterSolverAcc	IterSolver;


#else
#error "Incorrect user-supplied value for SOLVER. Check your build.config script."
#endif




#endif /* SOLVER_SPARSE_ITERSOLVERS_H_ */
