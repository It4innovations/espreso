
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPCG_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPCG_H_

#include "hypreboomeramg.h"
#include "hypreparasails.h"
#include "hypreeuclid.h"
#include "hyprepilut.h"

namespace espreso {

struct HYPREPCGConfiguration: public ECFDescription {

	enum class PRECONDITIONER {
		BoomerAMG,
		ParaSails,
		Euclid,
		Pilut,
		NONE
	};
	PRECONDITIONER preconditioner;

	HYPREBoomerAMGConfiguration boomeramg;
	HYPREParaSailsConfiguration parasails;
	HYPREEuclidConfiguration euclid;
	HYPREPilutConfiguration pilut;
	
	double relative_conv_tol, absolute_conv_tol, residual_conv_tol;
	int max_iterations;
	bool two_norm, recompute_residual_end, recompute_residual_p;
	
	enum class SOLVER_INFO {
		NO_INFO,
		SETUP_INFO,
		SOLVE_INFO,
		SETUP_SOLVE_INFO
	};
	SOLVER_INFO solver_info;

	HYPREPCGConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPCG_H_ */
