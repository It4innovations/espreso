
#ifndef SOLVER_SPECIFIC_DENSESOLVERS_H_
#define SOLVER_SPECIFIC_DENSESOLVERS_H_


#if defined(SOLVER_MKL)
#include "cpu/DenseSolverMKL.h"

namespace espreso {
	typedef DenseSolverMKL DenseSolverCPU;
	typedef DenseSolverMKL DenseSolverAcc;
}

#elif defined(SOLVER_PARDISO)
#include "cpu/solverpardiso.h"
#include "cpu/DenseSolverMKL.h"

namespace espreso {
	typedef DenseSolverMKL DenseSolverCPU;
	typedef DenseSolverMKL DenseSolverAcc;
}

#elif defined(SOLVER_MUMPS)
#include "cpu/solvermumps.h"

namespace espreso {
	typedef DenseSolverMUMPS DenseSolverCPU;
	typedef DenseSolverMUMPS DenseSolverAcc;
}

#elif defined(SOLVER_MIC)
#include "cpu/DenseSolverMKL.h"
//#include "acc/mic.h"

namespace espreso {
	typedef DenseSolverMKL DenseSolverCPU;
	typedef DenseSolverMKL DenseSolverAcc;
}

#elif defined(SOLVER_CUDA)
#include "cpu/DenseSolverMKL.h"
#include "acc/DenseSolverCUDA.h"

namespace espreso {
	typedef DenseSolverMKL DenseSolverCPU;
	typedef DenseSolverCUDA DenseSolverAcc;
}

#elif defined(SOLVER_CUDA_7)
#include "cpu/DenseSolverMKL.h"
#include "acc/DenseSolverCUDA.h"

namespace espreso {
	typedef DenseSolverMKL DenseSolverCPU;
	typedef DenseSolverCUDA DenseSolverAcc;
}


#else
#error "Incorrect user-supplied value for SOLVER. Check your build.config script."
#endif

#endif
