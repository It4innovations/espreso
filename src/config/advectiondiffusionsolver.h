
#ifndef SRC_CONFIG_ADVECTIONDIFFUSIONSOLVER_H_
#define SRC_CONFIG_ADVECTIONDIFFUSIONSOLVER_H_

#include "configuration.h"

namespace espreso {

struct NonLinearConvergence: public Configuration {

	PARAMETER(double, temperature, "Precision of temperature convergence.", 1e-3);
	PARAMETER(double, heat       , "Precision of heat convergence."       , 1e-3);

};

struct AutoTimeStepping: public Configuration {

	PARAMETER(bool, type, "Auto time stepping.", true);

	PARAMETER(double, initial_time_step, "Initial time step", 1e-2);
	PARAMETER(double, min_time_step, "Min time step", 1e-3);
	PARAMETER(double, max_time_step, "Max time step", 1);
};

struct AdvectionDiffusion2DNonLinearSolver: public Configuration {

	enum class METHOD {
		NEWTON_RHAPSON
	};

	OPTION(METHOD, method, "Non-linear method", METHOD::NEWTON_RHAPSON, OPTIONS({
		{ "NEWTON_RHAPSON", METHOD::NEWTON_RHAPSON, "Newton-Rhapson method." }
	}));

	PARAMETER(bool, line_search, "Set line search.", 1e-3);
	SUBCONFIG(NonLinearConvergence, convergence, "Convergence parameters.");
	PARAMETER(size_t, substeps, "Number of loading substeps.", 1);
};

struct AdvectionDiffusion2DTransientSolver: public Configuration {

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

struct AdvectionDiffusionSolver: public Configuration {

	enum class TYPE {
		STEADY_STATE,
		TRANSIENT
	};

	enum class MODE {
		LINEAR,
		NONLINEAR
	};

	enum class INTERPOLATION {
		LINEAR,
		QUADRATIC
	};

	OPTION(TYPE, type, "Solver type", TYPE::STEADY_STATE, OPTIONS({
		{ "STEADY_STATE", TYPE::STEADY_STATE, "Steady state." },
		{ "TRANSIENT"   , TYPE::TRANSIENT   , "Transient." },
	}));

	OPTION(MODE, mode, "Solver mode", MODE::LINEAR, OPTIONS({
		{ "LINEAR"   , MODE::LINEAR   , "Linear." },
		{ "NONLINEAR", MODE::NONLINEAR, "Non-linear." },
	}));

	OPTION(INTERPOLATION, interpolation, "Interpolation of tabular data.", INTERPOLATION::LINEAR, OPTIONS({
		{ "LINEAR"   , INTERPOLATION::LINEAR   , "Linear." },
		{ "QUADRATIC", INTERPOLATION::QUADRATIC, "Quadratic." },
	}));

	PARAMETER(size_t, load_steps, "Number of load steps in simulation.", 1);

	SUBMAPTOCONFIG(size_t, AdvectionDiffusion2DNonLinearSolver, nonlinear_solver, "Non-linear configuration for each load step.");
	SUBMAPTOCONFIG(size_t, AdvectionDiffusion2DTransientSolver, transient_solver, "Transient configuration for each load step.");
};

}


#endif /* SRC_CONFIG_ADVECTIONDIFFUSIONSOLVER_H_ */
