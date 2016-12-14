
#ifndef SRC_CONFIG_SOLVERHYPRE_H_
#define SRC_CONFIG_SOLVERHYPRE_H_

#include "configuration.h"

namespace espreso {

enum class HYPRE_SOLVER {
	CG = 0,
	GMRES = 1,
	FGMRES = 2,
	BICGS = 3,
	BICGSTAB = 4,
	TFQMR = 5,
	SYMQMR = 6,
	SUPERLU = 7,
	SUPERLUX = 8
};

enum class HYPRE_PRECONDITIONER {
	DIAGONAL = 0,
	PILUT = 1,
	EUCLID = 2,
	PARASAILS = 3,
	BOOMERAMG = 4,
	POLY = 5,
	MLI = 6
};

struct HypreSolver: public Configuration {

	PARAMETER(double, epsilon, "Solver requested precision.", 1e-5);
	PARAMETER(size_t, iterations, "solver max iterations", 100);

	OPTION(HYPRE_SOLVER, solver, "Used solver", HYPRE_SOLVER::CG, OPTIONS({
		{"CG"      , HYPRE_SOLVER::CG      , "CG solver." },
		{"GMRES"   , HYPRE_SOLVER::GMRES   , "GMRES solver." },
		{"FGMRES"  , HYPRE_SOLVER::FGMRES  , "FGMRES solver." },
		{"BICGS"   , HYPRE_SOLVER::BICGS   , "BICGS solver." },
		{"BICGSTAB", HYPRE_SOLVER::BICGSTAB, "BICGSTAB solver." },
		{"TFQMR"   , HYPRE_SOLVER::TFQMR   , "TFQMR solver." },
		{"SYMQMR"  , HYPRE_SOLVER::SYMQMR  , "SYMQMR solver." },
		{"SUPERLU" , HYPRE_SOLVER::SUPERLU , "SUPERLU solver." },
		{"SUPERLUX", HYPRE_SOLVER::SUPERLUX, "SUPERLUX solver." }
	}));

	OPTION(HYPRE_PRECONDITIONER, preconditioner, "Used preconditioner", HYPRE_PRECONDITIONER::BOOMERAMG, OPTIONS({
		{"DIAGONAL" , HYPRE_PRECONDITIONER::DIAGONAL , "DIAGONAL preconditioner." },
		{"PILUT"    , HYPRE_PRECONDITIONER::PILUT    , "PILUT preconditioner." },
		{"EUCLID"   , HYPRE_PRECONDITIONER::EUCLID   , "EUCLID preconditioner." },
		{"PARASAILS", HYPRE_PRECONDITIONER::PARASAILS, "PARASAILS preconditioner." },
		{"BOOMERAMG", HYPRE_PRECONDITIONER::BOOMERAMG, "BOOMERAMG preconditioner." },
		{"POLY"     , HYPRE_PRECONDITIONER::POLY     , "POLY preconditioner." },
		{"MLI"      , HYPRE_PRECONDITIONER::MLI      , "MLI preconditioner." }
	}));
};


}



#endif /* SRC_CONFIG_SOLVERHYPRE_H_ */
