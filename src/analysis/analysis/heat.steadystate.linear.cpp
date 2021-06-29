
#include "analysis.h"
#include "heat.steadystate.linear.h"
#include "analysis/linearsolver/directsolver.h"
#include "analysis/linearsolver/fetisolver.h"
#include "analysis/linearsolver/multigridsolver.h"

#include "config/ecf/physics/heattransfer.h"

using namespace espreso;

AX_HeatSteadyStateLinear::AX_HeatSteadyStateLinear(HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration)
: gsettings(gsettings), configuration(configuration), scheme{}, mode{}, solver{}
{
	switch (configuration.solver) {
	case LoadStepSolverConfiguration::SOLVER::FETI:    solver = Analysis::init<FETISolver<double> >(this, configuration.feti); break;
	case LoadStepSolverConfiguration::SOLVER::HYPRE:   solver = Analysis::init<MultigridSolver<double> >(this, configuration.hypre); break;
	case LoadStepSolverConfiguration::SOLVER::MKLPDSS: solver = Analysis::init<DirectSolver<double> >(this, configuration.mklpdss); break;
	case LoadStepSolverConfiguration::SOLVER::PARDISO: solver = Analysis::init<DirectSolver<double> >(this, configuration.pardiso); break;
	case LoadStepSolverConfiguration::SOLVER::SUPERLU: solver = Analysis::init<DirectSolver<double> >(this, configuration.superlu); break;
	case LoadStepSolverConfiguration::SOLVER::WSMP:    solver = Analysis::init<DirectSolver<double> >(this, configuration.wsmp); break;
	}
}

void AX_HeatSteadyStateLinear::solve()
{

}


