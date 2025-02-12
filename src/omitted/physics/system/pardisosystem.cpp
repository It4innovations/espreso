
#include "pardisosystem.h"
#include "builder/builder.h"

using namespace espreso;

PARDISOSystem::PARDISOSystem(int assemblers, int solvers, PARDISOConfiguration &configuration)
: assemblers(assemblers)
{
    this->solvers.reserve(solvers);
    for (int i = 0; i < solvers; i++) {
        this->solvers.emplace_back(configuration);
    }
}

void PARDISOSystem::_builderInit()
{
    builder->init(*this);
}

void PARDISOSystem::_builderReset()
{
    builder->reset(builder->matrices, *this);
}

void PARDISOSystem::_builderCreateSystem()
{
    builder->buildSystem(*this);
}

void PARDISOSystem::_builderUpdateSolution()
{
    builder->updateSolution(*this);
}


