
#ifndef SRC_CONFIGURATION_PHYSICS_TRANSIENTSOLVER_H_
#define SRC_CONFIGURATION_PHYSICS_TRANSIENTSOLVER_H_


#include "../configuration.h"

namespace espreso {

struct AutoTimeStepping: public Configuration {

	PARAMETER(bool, type, "Auto time stepping.", true);

	PARAMETER(double, initial_time_step, "Initial time step", 1e-2);
	PARAMETER(double, min_time_step, "Min time step", 1e-3);
	PARAMETER(double, max_time_step, "Max time step", 1);
};


struct TransientSolver: public Configuration {

	enum class METHOD {
		CRANK_NICOLSON,
		FORWARD_DIFF,
		GALERKIN,
		BACKWARD_DIFF
	};

	enum class ALPHA {
		OFF,
		ON,
		FROM_METHOD
	};

	OPTION(METHOD, method, "Transient method", METHOD::CRANK_NICOLSON, OPTIONS({
		{ "CRANK_NICOLSON", METHOD::CRANK_NICOLSON, "Crank-Nicolson method." },
		{ "FORWARD_DIFF"  , METHOD::FORWARD_DIFF  , "Forward differences." },
		{ "GALERKIN"      , METHOD::GALERKIN      , "Galerkin method." },
		{ "BACKWARD_DIFF" , METHOD::BACKWARD_DIFF , "Backward differences." }
	}));

	OPTION(ALPHA, alpha_set, "Alpha set", ALPHA::FROM_METHOD, OPTIONS({
		{ "0"          , ALPHA::OFF        , "Off." },
		{ "1"          , ALPHA::ON         , "On." },
		{ "FROM_METHOD", ALPHA::FROM_METHOD, "From method." }
	}));

	SUBCONFIG(AutoTimeStepping, auto_time_stepping, "Auto-time stepping parameters.");

	PARAMETER(double, time_step     , "Time step", 1e-2);
	PARAMETER(double, time_of_lading, "Lading time", 1);
};

}

#endif /* SRC_CONFIGURATION_PHYSICS_TRANSIENTSOLVER_H_ */
