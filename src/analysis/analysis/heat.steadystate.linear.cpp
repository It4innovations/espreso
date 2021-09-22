
#include "analysis.h"
#include "heat.steadystate.linear.h"

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

AX_HeatSteadyStateLinear::AX_HeatSteadyStateLinear(HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration)
: gsettings(gsettings), configuration(configuration), assembler{nullptr, gsettings, configuration}, scheme{}, system{}
{

}

AX_HeatSteadyStateLinear::~AX_HeatSteadyStateLinear()
{
	if (system) {
		delete system;
	}
}

void AX_HeatSteadyStateLinear::init()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info(" == ANALYSIS                                                            LINEAR STEADY STATE == \n");

//	switch (configuration.solver) {
//	case LoadStepSolverConfiguration::SOLVER::FETI:    system = new AX_FETISystem<double>(configuration.feti); break;
//	case LoadStepSolverConfiguration::SOLVER::HYPRE:   system = new AX_MultigridSystem<double>(configuration.hypre); break;
//	case LoadStepSolverConfiguration::SOLVER::MKLPDSS: system = new AX_MKLPDSSSystem<double>(configuration.mklpdss); break;
//	case LoadStepSolverConfiguration::SOLVER::PARDISO: system = new AX_DirectSystem<double>(configuration.pardiso); break;
//	case LoadStepSolverConfiguration::SOLVER::SUPERLU: system = new AX_DirectSystem<double>(configuration.superlu); break;
//	case LoadStepSolverConfiguration::SOLVER::WSMP:    system = new AX_DirectSystem<double>(configuration.wsmp); break;
//	}
	system = new AX_MKLPDSSSystem<AX_HeatSteadyStateLinear>(this, configuration.mklpdss);
	scheme.init(system);
	assembler.init(scheme);

	Variable::list.global.insert(std::make_pair("TIME", nullptr));

	info::mesh->output->updateMonitors(step::TYPE::TIME);
}

void AX_HeatSteadyStateLinear::run(step::Step &step)
{
	step::Time time;
	scheme.setTime(time, configuration.duration_time);
	Variable::list.global["TIME"] = new TimeVariable(time);

	assembler.evaluate();
	scheme.composeSystem(step, system);

	system->set(step);
	system->update(step);
	system->solve(step);

	scheme.extractSolution(step, system);

	assembler.updateSolution();
	info::mesh->output->updateSolution(step, time);

	eslog::info(" ============================================================================================= \n");
	eslog::info(" = ================================================================= run time %12.3f s =\n\n", eslog::duration());
}


