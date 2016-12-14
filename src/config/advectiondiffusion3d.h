
#ifndef SRC_CONFIG_ADVECTIONDIFFUSION3D_H_
#define SRC_CONFIG_ADVECTIONDIFFUSION3D_H_

#include "material.h"
#include "solver.h"

namespace espreso {

struct AdvectionDiffusion3DConfiguration: public Configuration {

	OPTION(SOLVER_LIBRARY, solver_library, "Linear solver used for computing a system.", SOLVER_LIBRARY::ESPRESO, OPTIONS({
		{ "ESPRESO", SOLVER_LIBRARY::ESPRESO, "ESPRESO solver [FETI methods]" },
		{ "HYPRE"  , SOLVER_LIBRARY::HYPRE  , "Hypre solver [multigrid methods]" },
	}));

	SUBCONFIG(ESPRESOSolver, espreso, "Internal FETI solver options.");
	SUBCONFIG(HypreSolver  , hypre  , "Multigrid solver setting.");
};

}



#endif /* SRC_CONFIG_ADVECTIONDIFFUSION3D_H_ */
