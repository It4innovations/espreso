
#include "analysis.h"
#include "heat.steadystate.nonlinear.h"

#include "analysis/linearsystem/directsystem.h"
#include "analysis/linearsystem/fetisystem.h"
#include "analysis/linearsystem/mklpdsssystem.h"
#include "analysis/linearsystem/multigridsystem.h"

#include "basis/expression/variable.h"
#include "config/ecf/physics/heattransfer.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "output/output.h"

using namespace espreso;

AX_HeatSteadyStateNonLinear::AX_HeatSteadyStateNonLinear(HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration)
: gsettings(gsettings), configuration(configuration), assembler{nullptr, gsettings, configuration}, solver(configuration.nonlinear_solver), scheme{}, system{}

{

}

AX_HeatSteadyStateNonLinear::~AX_HeatSteadyStateNonLinear()
{
	if (system) {
		delete system;
	}
}

void AX_HeatSteadyStateNonLinear::init()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info(" == ANALYSIS                                                        NON-LINEAR STEADY STATE == \n");

//	switch (configuration.solver) {
//	case LoadStepSolverConfiguration::SOLVER::FETI:    system = new AX_FETISystem<double>(configuration.feti); break;
//	case LoadStepSolverConfiguration::SOLVER::HYPRE:   system = new AX_MultigridSystem<double>(configuration.hypre); break;
//	case LoadStepSolverConfiguration::SOLVER::MKLPDSS: system = new AX_MKLPDSSSystem<double>(configuration.mklpdss); break;
//	case LoadStepSolverConfiguration::SOLVER::PARDISO: system = new AX_DirectSystem<double>(configuration.pardiso); break;
//	case LoadStepSolverConfiguration::SOLVER::SUPERLU: system = new AX_DirectSystem<double>(configuration.superlu); break;
//	case LoadStepSolverConfiguration::SOLVER::WSMP:    system = new AX_DirectSystem<double>(configuration.wsmp); break;
//	}
	system = new AX_MKLPDSSSystem<AX_HeatSteadyStateNonLinear>(this, configuration.mklpdss);
	solver.init(system);
	scheme.init(system);
	assembler.init(scheme, system->assembler.dirichlet);

	Variable::list.global.insert(std::make_pair("TIME", nullptr));

	info::mesh->output->updateMonitors(step::TYPE::TIME);
}

void AX_HeatSteadyStateNonLinear::run(step::Step &step)
{
	step::Time time;
	scheme.setTime(time, configuration.duration_time);
	Variable::list.global["TIME"] = new TimeVariable(time);

	solver.run(step, time, assembler, scheme, system);

	info::mesh->output->updateSolution(time);

	eslog::info(" ============================================================================================= \n");
	eslog::info(" =================================================================== run time %12.3f s =\n\n", eslog::duration());
}


