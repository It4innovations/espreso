
#include "superlusystem.h"
#include "builder/builder.h"

using namespace espreso;

SuperLUSystem::SuperLUSystem(int assemblers, int solvers, SuperLUConfiguration &configuration)
: assemblers(assemblers)
{
    this->solvers.reserve(solvers);
    for (int i = 0; i < solvers; i++) {
        this->solvers.emplace_back(configuration);
    }
}

void SuperLUSystem::_builderInit()
{
    builder->init(*this);
}

void SuperLUSystem::_builderReset()
{
    builder->reset(builder->matrices, *this);
}

void SuperLUSystem::_builderCreateSystem()
{
    builder->buildSystem(*this);
}

void SuperLUSystem::_builderUpdateSolution()
{
    builder->updateSolution(*this);
}
