
#ifndef SOLVER_SPECIFIC_ITERSOLVERS_H_
#define SOLVER_SPECIFIC_ITERSOLVERS_H_

#if defined(SOLVER_MKL)
#include "cpu/itersolvercpu.h"

namespace espreso {
	typedef IterSolverCPU	IterSolver;
}


#elif defined(SOLVER_PARDISO)
#include "cpu/itersolvercpu.h"

namespace espreso {
	typedef IterSolverCPU	IterSolver;
}

#elif defined(SOLVER_MUMPS)
#include "cpu/itersolvercpu.h"

namespace espreso {
	typedef IterSolverCPU	IterSolver;
}

#elif defined(SOLVER_MIC)
#include "acc/itersolveracc.h"

namespace espreso {
	typedef IterSolverAcc	IterSolver;
}

#elif defined(SOLVER_CUDA)
#include "acc/itersolverGPU.h"

namespace espreso {
	typedef IterSolverGPU	IterSolver;
}

#elif defined(SOLVER_CUDA_7)
#include "acc/itersolverGPU.h"

namespace espreso {
	typedef IterSolverGPU	IterSolver;
}

#elif defined(SOLVER_DISSECTION)
#include "cpu/itersolvercpu.h"

namespace espreso {
	typedef IterSolverCPU	IterSolver;
}



#else
#error "Incorrect user-supplied value for SOLVER. Check your build.config script."
#endif




#endif /* SOLVER_SPECIFIC_ITERSOLVERS_H_ */
