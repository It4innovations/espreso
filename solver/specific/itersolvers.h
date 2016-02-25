
#ifndef SOLVER_SPECIFIC_ITERSOLVERS_H_
#define SOLVER_SPECIFIC_ITERSOLVERS_H_

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
#include "acc/itersolverGPU.h"
	typedef IterSolverGPU	IterSolver;


#else
#error "Incorrect user-supplied value for SOLVER. Check your build.config script."
#endif




#endif /* SOLVER_SPECIFIC_ITERSOLVERS_H_ */
