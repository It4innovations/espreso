
#ifndef SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_
#define SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_

#include "nonlinear.h"
#include "transientfirstorderimplicit.h"
#include "transientsecondorderimplicit.h"
#include "harmonic.h"
#include "config/ecf/linearsolver/feti.h"
#include "config/ecf/linearsolver/hypre/hypre.h"
#include "config/ecf/linearsolver/mklpdss.h"
#include "config/ecf/linearsolver/pardiso.h"
#include "config/ecf/linearsolver/superlu.h"
#include "config/ecf/linearsolver/wsmp.h"

namespace espreso {

struct TopologyOptimizationOCSolverSettings: public ECFDescription {
	double lower_bound;
	double upper_bound;
	double move;
	double precision;

	TopologyOptimizationOCSolverSettings();
};

struct TopologyOptimizationMMASolverSettings: public ECFDescription {
	TopologyOptimizationMMASolverSettings();
};

struct TopologyOptimizationSolverSettings: public ECFDescription {
	enum class Type {
		OC,
		MMA
	};

	Type type;

	double min_density;
	double precision;
	int max_iterations;
	double penalty_factor;

	TopologyOptimizationOCSolverSettings oc;
	TopologyOptimizationMMASolverSettings mma;

	TopologyOptimizationSolverSettings();
};

struct TopologyOptimizationConstraint: public ECFDescription {
	enum class Response {
		VOLUME,
		MASS
	};

	enum class PresetValues {
		SOLID,
		VOID
	};

	Response response;
	double value;
	std::map<std::string, PresetValues> preset_regions;

	TopologyOptimizationConstraint();
};

struct TopologyOptimizationFilteringDensity: public ECFDescription {
	enum class Type {
		LINEAR,
		GAUSSIAN,
		HEAVISIDE
	};

	Type type;

	TopologyOptimizationFilteringDensity();
};

struct TopologyOptimizationFiltering: public ECFDescription {
	double radius;
	TopologyOptimizationFilteringDensity density;

	TopologyOptimizationFiltering();
};

struct TopologyOptimizationConfiguration: public ECFDescription {

	TopologyOptimizationSolverSettings solver_settings;
	TopologyOptimizationConstraint constraint;
	TopologyOptimizationFiltering filtering;

	TopologyOptimizationConfiguration();
};

struct LoadStepSolverConfiguration: public ECFDescription {

	enum class TYPE {
		STEADY_STATE,
		TRANSIENT,
		HARMONIC
	};

	enum class MODE {
		LINEAR,
		NONLINEAR
	};

	enum class SOLVER {
		FETI,
		HYPRE,
		MKLPDSS,
		PARDISO,
		SUPERLU,
		WSMP
	};

	double duration_time;
	bool topology_optimization;

	TYPE type;
	MODE mode;
	SOLVER solver;

	TopologyOptimizationConfiguration topology_optimization_settings;

	FETIConfiguration feti;
	HYPREConfiguration hypre;
	MKLPDSSConfiguration mklpdss;
	PARDISOConfiguration pardiso;
	SuperLUConfiguration superlu;
	WSMPConfiguration wsmp;

	LoadStepSolverConfiguration();
};

struct HeatTransferLoadStepSolverConfiguration: public LoadStepSolverConfiguration {

	NonLinearSolverConfiguration nonlinear_solver;
	TransientFirstOrderImplicitSolverConfiguration transient_solver;

	HeatTransferLoadStepSolverConfiguration();
};

struct AcousticLoadStepSolverConfiguration: public LoadStepSolverConfiguration {

	HarmonicSolverConfiguration harmonic_solver;

	AcousticLoadStepSolverConfiguration();
};

struct StructuralMechanicsLoadStepSolverConfiguration: public LoadStepSolverConfiguration {

	NonLinearSolverConfiguration nonlinear_solver;
	TransientSecondOrderImplicitSolverConfiguration transient_solver;
	HarmonicSolverConfiguration harmonic_solver;

	StructuralMechanicsLoadStepSolverConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_ */
