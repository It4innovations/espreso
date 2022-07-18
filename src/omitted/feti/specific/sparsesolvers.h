
#ifndef SOLVER_SPECIFIC_SPARSESOLVERS_H_
#define SOLVER_SPECIFIC_SPARSESOLVERS_H_


#if defined(SOLVER_MKL)
#include "cpu/SparseSolverMKL.h"

namespace espreso {
	typedef SparseSolverMKL SparseSolverCPU;
	typedef SparseSolverMKL SparseSolverAcc;
}


#elif defined(SOLVER_PARDISO)
#include "cpu/SparseSolverPARDISO.h"

namespace espreso {
	typedef SparseSolverPARDISO SparseSolverCPU;
	typedef SparseSolverPARDISO SparseSolverAcc;
}

#elif defined(SOLVER_MUMPS)
#include "cpu/solvermumps.h"

namespace espreso {
	typedef SparseSolverMUMPS SparseSolverCPU;
	typedef SparseSolverMUMPS SparseSolverAcc;
}

#elif defined(SOLVER_MIC)
#include "cpu/SparseSolverMKL.h"
#include "acc/mic.h"

namespace espreso {
	typedef SparseSolverMKL SparseSolverCPU;
	typedef SparseSolverMIC SparseSolverAcc;
}


#elif defined(SOLVER_CUDA_X)
#include "cpu/SparseSolverMKL.h"
#include "acc/SparseSolverCUDA.h"

namespace espreso {
	typedef SparseSolverMKL SparseSolverCPU;
	typedef SparseSolverCUDA SparseSolverAcc;
}

#elif defined(SOLVER_CUDA)
#include "cpu/SparseSolverMKL.h"

namespace espreso {
	typedef SparseSolverMKL SparseSolverCPU;
	typedef SparseSolverMKL SparseSolverAcc;
}

#elif defined(SOLVER_DISSECTION)
#include "cpu/SparseSolverDissection.h"
#include "cpu/SparseSolverMKL.h"

namespace espreso {
	typedef SparseSolverDissection SparseSolverCPU;
	typedef SparseSolverDissection SparseSolverAcc;
}


#else
#error "Incorrect user-supplied value for SOLVER. Check your build.config script."
#endif




#endif /* SOLVER_SPECIFIC_SPARSESOLVERS_H_ */
