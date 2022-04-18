
#include "mklpdsssystem.h"
#include "builder/builder.h"

using namespace espreso;

MKLPDSSSystem::MKLPDSSSystem(int assemblers, int solvers, MKLPDSSConfiguration &configuration)
: assemblers(assemblers)
{
	this->solvers.reserve(solvers);
	for (int i = 0; i < solvers; i++) {
		this->solvers.emplace_back(configuration);
	}
}

void MKLPDSSSystem::_builderInit()
{
	builder->init(*this);
}

void MKLPDSSSystem::_builderReset()
{
	builder->reset(builder->matrices, *this);
}

void MKLPDSSSystem::_builderCreateSystem()
{
	builder->buildSystem(*this);
}

void MKLPDSSSystem::_builderUpdateSolution()
{
	builder->updateSolution(*this);
}
