
#include "analysis.h"
#include "heat.steadystate.linear.h"

#include "analysis/linearsystem/directsystem.h"
#include "analysis/linearsystem/fetisystem.h"
#include "analysis/linearsystem/mklpdsssystem.h"
#include "analysis/linearsystem/multigridsystem.h"

#include "config/ecf/physics/heattransfer.h"

#include "esinfo/meshinfo.h"
#include "output/output.h"

using namespace espreso;

AX_HeatSteadyStateLinear::AX_HeatSteadyStateLinear(HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration)
: gsettings(gsettings), configuration(configuration), assembler{nullptr, gsettings, configuration}, scheme{}, system{}
{

}

void AX_HeatSteadyStateLinear::init()
{
	switch (configuration.solver) {
	case LoadStepSolverConfiguration::SOLVER::FETI:    system = new AX_FETISystem<double>(configuration.feti); break;
	case LoadStepSolverConfiguration::SOLVER::HYPRE:   system = new AX_MultigridSystem<double>(configuration.hypre); break;
	case LoadStepSolverConfiguration::SOLVER::MKLPDSS: system = new AX_MKLPDSSSystem<double>(configuration.mklpdss); break;
	case LoadStepSolverConfiguration::SOLVER::PARDISO: system = new AX_DirectSystem<double>(configuration.pardiso); break;
	case LoadStepSolverConfiguration::SOLVER::SUPERLU: system = new AX_DirectSystem<double>(configuration.superlu); break;
	case LoadStepSolverConfiguration::SOLVER::WSMP:    system = new AX_DirectSystem<double>(configuration.wsmp); break;
	}
	system->init(this);
	scheme.init(system);
	assembler.init(scheme);

	info::mesh->output->updateMonitors(step::TYPE::TIME);
}

void AX_HeatSteadyStateLinear::run()
{
	step::Time time;
	scheme.setTime(time, configuration.duration_time);

	assembler.next();
	scheme.composeSystem(system);
	assembler.fillDirichlet(*system->solver.dirichlet);

	scheme.storeScheme(time);

	system->update(assembler);
	system->solve();

	scheme.extractSolution(system);
	scheme.storeSolution(time);

	assembler.updateSolution();
	info::mesh->output->updateSolution(time);
}


