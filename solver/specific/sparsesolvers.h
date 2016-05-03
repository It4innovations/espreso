
#ifndef SOLVER_SPECIFIC_SPARSESOLVERS_H_
#define SOLVER_SPECIFIC_SPARSESOLVERS_H_


#if defined(SOLVER_MKL)
#include "cpu/solvermkl.h"

namespace espreso {
	typedef SparseSolverMKL SparseSolverCPU;
	typedef SparseSolverMKL SparseSolverAcc;
}


#elif defined(SOLVER_PARDISO)
#include "cpu/solverpardiso.h"

namespace espreso {
	typedef SparseSolverPardiso SparseSolverCPU;
	typedef SparseSolverPardiso SparseSolverAcc;
}

#elif defined(SOLVER_MUMPS)
#include "cpu/solvermumps.h"

namespace espreso {
	typedef SparseSolverMUMPS SparseSolverCPU;
	typedef SparseSolverMUMPS SparseSolverAcc;
}

#elif defined(SOLVER_MIC)
#include "cpu/solvermkl.h"
#include "acc/mic.h"

namespace espreso {
	typedef SparseSolverMKL SparseSolverCPU;
	typedef SparseSolverMIC SparseSolverAcc;
}


#elif defined(SOLVER_CUDA)
#include "cpu/solvermkl.h"
//#include "acc/solvercuda.h"

namespace espreso {
	typedef SparseSolverMKL SparseSolverCPU;
	typedef SparseSolverMKL SparseSolverAcc;
	//typedef SparseSolverCUDA SparseSolverAcc;
}


#else
#error "Incorrect user-supplied value for SOLVER. Check your build.config script."
#endif




#endif /* SOLVER_SPECIFIC_SPARSESOLVERS_H_ */
