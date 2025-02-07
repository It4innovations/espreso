
#include "loadstep.h"
#include "config/configuration.hpp"

using namespace espreso;

TopologyOptimizationOCSolverSettings::TopologyOptimizationOCSolverSettings()
{
	lower_bound = 0;
	REGISTER(lower_bound, ECFMetaData()
			.setdescription({ "Lower bounds" })
			.setdatatype({ ECFDataType::FLOAT }));
	upper_bound = 1e9;
	REGISTER(upper_bound, ECFMetaData()
			.setdescription({ "Upper bounds" })
			.setdatatype({ ECFDataType::FLOAT }));
	move = .2;
	REGISTER(move, ECFMetaData()
			.setdescription({ "Move" })
			.setdatatype({ ECFDataType::FLOAT }));
	precision = 1e-3;
	REGISTER(precision, ECFMetaData()
			.setdescription({ "Precision" })
			.setdatatype({ ECFDataType::FLOAT }));
}

TopologyOptimizationMMASolverSettings::TopologyOptimizationMMASolverSettings()
{

}

TopologyOptimizationSolverSettings::TopologyOptimizationSolverSettings()
{
	type = Type::OC;
	REGISTER(type, ECFMetaData()
			.setdescription({ "Solver type." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("OC").setdescription("..."))
			.addoption(ECFOption().setname("MMA").setdescription("...")));

	min_density = 1e-5;
	REGISTER(min_density, ECFMetaData()
			.setdescription({ "Min density" })
			.setdatatype({ ECFDataType::FLOAT }));
	precision = 1e-2;
	REGISTER(precision, ECFMetaData()
			.setdescription({ "Precision" })
			.setdatatype({ ECFDataType::FLOAT }));
	max_iterations = 100;
	REGISTER(max_iterations, ECFMetaData()
			.setdescription({ "Max iterations." })
			.setdatatype({ ECFDataType::INTEGER }));
	penalty_factor = 3;
	REGISTER(penalty_factor, ECFMetaData()
			.setdescription({ "Penalty factor." })
			.setdatatype({ ECFDataType::FLOAT }));

	REGISTER(oc, ECFMetaData()
		.setdescription({ "OC solver settings." })
		.setcollapsed());
	REGISTER(mma, ECFMetaData()
		.setdescription({ "MMA solver settings." })
		.setcollapsed());
}

TopologyOptimizationConstraint::TopologyOptimizationConstraint()
{
	response = Response::VOLUME;
	REGISTER(response, ECFMetaData()
			.setdescription({ "Response type." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("VOLUME").setdescription("Steady state load step."))
			.addoption(ECFOption().setname("MASS").setdescription("Transient load step.")));

	value = .5;
	REGISTER(value, ECFMetaData()
			.setdescription({ "Value." })
			.setdatatype({ ECFDataType::FLOAT }));

	REGISTER(preset_regions, ECFMetaData()
		.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::OPTION })
		.setdescription({ "Regions with fixed design variable.", "Preset value." })
		.setdynamic()
		.setpattern({ "MY_REGION", "SOLID" })
		.addoption(ECFOption().setname("SOLID").setdescription("Solid region."))
		.addoption(ECFOption().setname("VOID").setdescription("Void region.")));
}

TopologyOptimizationFilteringDensity::TopologyOptimizationFilteringDensity()
{
	type = Type::LINEAR;
	REGISTER(type, ECFMetaData()
			.setdescription({ "Density type." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("LINEAR").setdescription("Linear."))
			.addoption(ECFOption().setname("GAUSSIAN").setdescription("Gaussian."))
			.addoption(ECFOption().setname("HEAVISIDE").setdescription("Heaviside.")));
}

TopologyOptimizationFiltering::TopologyOptimizationFiltering()
{
	radius = 1;
	REGISTER(radius, ECFMetaData()
			.setdescription({ "Radius." })
			.setdatatype({ ECFDataType::FLOAT }));
	REGISTER(density, ECFMetaData()
		.setdescription({ "Filtering density." })
		.setgroup());
}

TopologyOptimizationConfiguration::TopologyOptimizationConfiguration()
{
	REGISTER(solver_settings, ECFMetaData()
		.setdescription({ "Topology optimization solver settings." })
		.setcollapsed());
	
	REGISTER(constraint, ECFMetaData()
		.setdescription({ "Topology optimization constraint." })
		.setcollapsed());
	
	REGISTER(filtering, ECFMetaData()
		.setdescription({ "Topology optimization filtering." })
		.setcollapsed());
}

LoadStepSolverConfiguration::LoadStepSolverConfiguration()
{
	duration_time = 1;
	REGISTER(duration_time, ECFMetaData()
			.setdescription({ "Duration" })
			.setdatatype({ ECFDataType::FLOAT }));

	topology_optimization = false;
	REGISTER(topology_optimization, ECFMetaData()
			.setdescription({ "Turn ON/OFF topology optimization." })
			.setdatatype({ ECFDataType::BOOL }));

	type = TYPE::STEADY_STATE;
	REGISTER(type, ECFMetaData()
			.setdescription({ "Simulation type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("STEADY_STATE").setdescription("Steady state load step."))
			.addoption(ECFOption().setname("TRANSIENT").setdescription("Transient load step."))
			.addoption(ECFOption().setname("HARMONIC").setdescription("Harmonic load step.")));

	mode = MODE::LINEAR;
	REGISTER(mode, ECFMetaData()
			.setdescription({ "Simulation mode" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("LINEAR").setdescription("Linear material behavior."))
			.addoption(ECFOption().setname("NONLINEAR").setdescription("Nonlinear material behavior.")));

	solver = SOLVER::FETI;
	REGISTER(solver, ECFMetaData()
			.setdescription({ "Linear solver" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("FETI").setdescription("Use ESPRESO as linear solver."))
			.addoption(ECFOption().setname("HYPRE").setdescription("Use hypre library."))
			.addoption(ECFOption().setname("MKLPDSS").setdescription("Use parallel direct sparse solver from MKL (PARDISO)."))
			.addoption(ECFOption().setname("PARDISO").setdescription("Use PARDISO original library."))
			.addoption(ECFOption().setname("SUPERLU").setdescription("Use SuperLU_DIST."))
			.addoption(ECFOption().setname("SEQUENTIAL").setdescription("Use sequential SpSolver."))
			.addoption(ECFOption().setname("WSMP").setdescription("Use Watson Sparse Matrix Package."))
			.addoption(ECFOption().setname("NONE").setdescription("Dry run just with assembler.")));

	REGISTER(topology_optimization_settings, ECFMetaData()
			.setdescription({ "Topology optimization settings." })
			.setcollapsed()
			.allowonly([&] () { return topology_optimization; }));

	REGISTER(feti, ECFMetaData()
			.setdescription({ "FETI solver settings" })
			.setcollapsed()
			.allowonly([&] () { return solver == SOLVER::FETI; }));
	REGISTER(hypre, ECFMetaData()
			.setdescription({ "HYPRE multigrid solver settings" })
			.setcollapsed()
			.allowonly([&] () { return solver == SOLVER::HYPRE; }));
	REGISTER(mklpdss, ECFMetaData()
			.setdescription({ "MKL parallel direct sparse solver" })
			.noexport()
			.allowonly([&] () { return solver == SOLVER::MKLPDSS; }));
	REGISTER(pardiso, ECFMetaData()
			.setdescription({ "PARDISO direct solver" })
			.noexport()
			.allowonly([&] () { return solver == SOLVER::PARDISO; }));
	REGISTER(superlu, ECFMetaData()
			.setdescription({ "SuperLU_DIST direct sparse solver" })
			.setcollapsed()
			.allowonly([&] () { return solver == SOLVER::SUPERLU; }));
	REGISTER(wsmp, ECFMetaData()
			.setdescription({ "WSMP direct sparse solver" })
			.noexport()
			.allowonly([&] () { return solver == SOLVER::WSMP; }));
}

HeatTransferLoadStepSolverConfiguration::HeatTransferLoadStepSolverConfiguration()
: nonlinear_solver("temperature", "heat")
{
	REGISTER(nonlinear_solver, ECFMetaData()
			.setdescription({ "Non-linear physics solver settings" })
			.setgroup()
			.allowonly([&] () { return mode == MODE::NONLINEAR; }));

	REGISTER(transient_solver, ECFMetaData()
			.setdescription({ "Transient physics solver settings" })
			.setgroup()
			.allowonly([&] () { return type == TYPE::TRANSIENT; }));
}

AcousticLoadStepSolverConfiguration::AcousticLoadStepSolverConfiguration()
{
	REGISTER(harmonic_solver, ECFMetaData()
			.setdescription({ "Harmonic physics solver settings" })
			.allowonly([&] () { return type == TYPE::HARMONIC; }));
}

StructuralMechanicsLoadStepSolverConfiguration::StructuralMechanicsLoadStepSolverConfiguration()
: nonlinear_solver("displacement", "forces")
{
	REGISTER(nonlinear_solver, ECFMetaData()
			.setdescription({ "Non-linear physics solver settings" })
			.allowonly([&] () { return mode == MODE::NONLINEAR; }));

	REGISTER(transient_solver, ECFMetaData()
			.setdescription({ "Transient physics solver settings" })
			.allowonly([&] () { return type == TYPE::TRANSIENT; }));

	REGISTER(harmonic_solver, ECFMetaData()
			.setdescription({ "Harmonic physics solver settings" })
			.allowonly([&] () { return type == TYPE::HARMONIC; }));
}


