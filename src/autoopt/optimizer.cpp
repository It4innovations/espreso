
#include "optimizer.h"
#include "esinfo/eslog.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"

#include <iostream>

using namespace espreso;

bool EmptyOptimizer::call(std::function<bool(void)> fnc)
{
	bool ret = fnc();
	if (!ret) { eslog::error("Called FETI with invalid parameters!\n"); }
	return ret;
}

EvolutionaryOptimizer::EvolutionaryOptimizer(const AutoOptimizationConfiguration& configuration, std::vector<ECFParameter*>& parameters)
: m_proxy(parameters, configuration)
{
	// sphere.forEachParameters(
	// 	[&] (ECFParameter* p) { this->addParameter(p); }
	// );
}


bool EvolutionaryOptimizer::call(std::function<bool(void)> fnc)
{
	this->m_proxy.setNextConfiguration();
	
	double start = eslog::time();
	int l_ret = fnc();
	double end = eslog::time();

	int g_ret = 0;
	Communication::allReduce(&l_ret, &g_ret, 1, MPI_INT, MPI_MIN);
	
	if (!g_ret)
	{
		this->m_proxy.setConfigurationForbidden();
	}
	else 
	{
		this->m_proxy.setConfigurationEvaluation(end - start);
	}

	return static_cast<bool>(g_ret);
	// for (auto p = _parameters.begin(); p != _parameters.end(); ++p) {
	// 	std::cout << (*p)->name << ": " << (*p)->getValue() << " ";
	// }
	// std::cout << std::endl;
}

