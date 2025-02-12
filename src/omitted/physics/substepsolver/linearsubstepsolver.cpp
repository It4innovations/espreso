
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
    _system->builder->matrices = Builder::Request::KCM | Builder::Request::f | Builder::Request::BC;
    loadStepSolver.updateStructuralMatrices();
    _system->setDirichlet();

    eslog::solver("     - LINEAR TIME STEP                          REASSEMBLED MATRICES  ::  %c, %c, %c, %c, %c -    \n",
            (_system->builder->matrices & Builder::Request::K)  ? 'K' : ' ',
            (_system->builder->matrices & Builder::Request::M)  ? 'M' : ' ',
            (_system->builder->matrices & Builder::Request::C)  ? 'C' : ' ',
            (_system->builder->matrices & Builder::Request::R)  ? 'R' : ' ',
            (_system->builder->matrices & Builder::Request::f)  ? 'f' : ' ');

    _system->solve();
    _system->solutionChanged();
}




