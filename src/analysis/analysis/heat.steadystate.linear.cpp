
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
: gsettings(gsettings), configuration(configuration), assembler{nullptr, gsettings, configuration}, scheme{}, mode{}, system{}
{
	switch (configuration.solver) {
	case LoadStepSolverConfiguration::SOLVER::FETI:    system = Analysis::init<AX_FETISystem, double>(this, configuration.feti); break;
	case LoadStepSolverConfiguration::SOLVER::HYPRE:   system = Analysis::init<AX_MultigridSystem, double>(this, configuration.hypre); break;
	case LoadStepSolverConfiguration::SOLVER::MKLPDSS: system = Analysis::init<AX_MKLPDSSSystem, double>(this, configuration.mklpdss); break;
	case LoadStepSolverConfiguration::SOLVER::PARDISO: system = Analysis::init<AX_DirectSystem, double>(this, configuration.pardiso); break;
	case LoadStepSolverConfiguration::SOLVER::SUPERLU: system = Analysis::init<AX_DirectSystem, double>(this, configuration.superlu); break;
	case LoadStepSolverConfiguration::SOLVER::WSMP:    system = Analysis::init<AX_DirectSystem, double>(this, configuration.wsmp); break;
	}
	info::mesh->output->updateMonitors();
}

void AX_HeatSteadyStateLinear::solve()
{
	if (mode.solve(scheme, assembler)) {
		scheme.solved();
	} else {
		// error
	}
}


