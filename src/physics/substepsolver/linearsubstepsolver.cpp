
#include "linearsubstepsolver.h"
#include "esinfo/eslog.hpp"
#include "physics/system/linearsystem.h"
#include "physics/system/builder/builder.h"
#include "physics/loadstepsolver/loadstepsolver.h"

#include "config/ecf/physics/physicssolver/loadstep.h"

using namespace espreso;

LinearSubStep::LinearSubStep(LinearSystem *system)
: SubStepSolver(system)
{

}

void LinearSubStep::init(SubStepSolver *previous)
{

}

bool LinearSubStep::hasSameMode(const LoadStepSolverConfiguration &configuration) const
{
	return configuration.mode == LoadStepSolverConfiguration::MODE::LINEAR;
}

void LinearSubStep::solve(LoadStepSolver &loadStepSolver)
{
	system->builder->matrices = Builder::Request::KCM | Builder::Request::f | Builder::Request::BC;
	loadStepSolver.updateStructuralMatrices();
	system->setDirichlet();

	eslog::solver("     - LINEAR TIME STEP                          REASSEMBLED MATRICES  ::  %c, %c, %c, %c, %c -    \n",
			(system->builder->matrices & Builder::Request::K)  ? 'K' : ' ',
			(system->builder->matrices & Builder::Request::M)  ? 'M' : ' ',
			(system->builder->matrices & Builder::Request::C)  ? 'C' : ' ',
			(system->builder->matrices & Builder::Request::R)  ? 'R' : ' ',
			(system->builder->matrices & Builder::Request::f)  ? 'f' : ' ');

	system->solve();
	system->solutionChanged();
}




