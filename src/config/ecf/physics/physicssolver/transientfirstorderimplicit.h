
#ifndef SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_TRANSIENTFIRSTORDERIMPLICIT_H_
#define SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_TRANSIENTFIRSTORDERIMPLICIT_H_

#include "config/description.h"

namespace espreso {

struct AutoTimeSteppingConfiguration: public ECFDescription {

	bool allowed;

	double min_time_step, max_time_step;

	double oscilation_limit, IDFactor;
	int points_per_period;

	AutoTimeSteppingConfiguration();
};


struct TransientFirstOrderImplicitSolverConfiguration: public ECFDescription {

	enum class METHOD {
		CRANK_NICOLSON,
		FORWARD_DIFF,
		GALERKIN,
		BACKWARD_DIFF,
		USER
	};

	METHOD method;
	AutoTimeSteppingConfiguration auto_time_stepping;
	double alpha, time_step;

	TransientFirstOrderImplicitSolverConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_TRANSIENTFIRSTORDERIMPLICIT_H_ */
