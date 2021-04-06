
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPRECGNR_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPRECGNR_H_

#include "hypreboomeramg.h"

namespace espreso {

struct HYPRECGNRConfiguration: public ECFDescription {

	enum class PRECONDITIONER {
		BoomerAMG,
	
		NONE
	};
	PRECONDITIONER preconditioner;

	HYPREBoomerAMGConfiguration boomeramg;

	double relative_conv_tol, absolute_conv_tol;
	int max_iterations;
	
	enum class SOLVER_INFO {
		NO_INFO,
		SETUP_INFO,
		SOLVE_INFO,
		SETUP_SOLVE_INFO
	};
	SOLVER_INFO solver_info;

	HYPRECGNRConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPRECGNR_H_ */
