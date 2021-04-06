
#include "hypresystem.h"
#include "builder/builder.h"
#include "basis/utilities/debugprint.h"
#include "basis/utilities/sysutils.h"

#include <fstream>

using namespace espreso;

void HYPREAssemblerData::print(const Builder *builder, const char* prefix, const char* suffix)
{
	DistributedAssemblerData::print(builder, prefix, suffix);
	if (builder->matrices & Builder::Request::K) {
		std::ofstream os(utils::prepareFile(std::string(prefix), std::string("N") + std::string(suffix)));
		os << N;
	}
}

void HYPRESolverData::printData(const Builder *builder, const char* prefix)
{
	DistributedSolverData::printData(builder, prefix);
	if (builder->matrices & Builder::Request::K) {
		std::ofstream os(utils::prepareFile(std::string(prefix), std::string("N")));
		os << N;
	}
}

HYPRESystem::HYPRESystem(int assemblers, int solvers, HYPREConfiguration &configuration)
: assemblers(assemblers)
{
	this->solvers.reserve(solvers);
	for (int i = 0; i < solvers; i++) {
		this->solvers.emplace_back(configuration);
	}
}

void HYPRESystem::_builderInit()
{
	builder->init(*this);
}

void HYPRESystem::_builderReset()
{
	builder->reset(builder->matrices, *this);
}

void HYPRESystem::_builderCreateSystem()
{
	builder->buildSystem(*this);
}

void HYPRESystem::_builderUpdateSolution()
{
	builder->updateSolution(*this);
}
