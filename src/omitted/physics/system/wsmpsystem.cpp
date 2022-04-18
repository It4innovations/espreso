
#include "wsmpsystem.h"
#include "builder/builder.h"

using namespace espreso;

WSMPSystem::WSMPSystem(int assemblers, int solvers, WSMPConfiguration &configuration)
: assemblers(assemblers)
{
	this->solvers.reserve(solvers);
	for (int i = 0; i < solvers; i++) {
		this->solvers.emplace_back(configuration);
	}
}

void WSMPSystem::_builderInit()
{
	builder->init(*this);
}

void WSMPSystem::_builderReset()
{
	builder->reset(builder->matrices, *this);
}

void WSMPSystem::_builderCreateSystem()
{
	builder->buildSystem(*this);
}

void WSMPSystem::_builderUpdateSolution()
{
	builder->updateSolution(*this);
}
