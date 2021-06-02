
#include "optimizer.h"
#include "esinfo/eslog.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"

#include <iostream>

using namespace espreso;

bool EmptyOptimizer::set(std::function<bool(void)> fnc)
{
	return fnc();
}

bool EmptyOptimizer::run(std::function<bool(void)> fnc)
{
	return fnc();
}

EvolutionaryOptimizer::EvolutionaryOptimizer(const AutoOptimizationConfiguration& configuration, std::vector<ECFParameter*>& parameters)
: m_proxy(parameters, configuration)
{
	// sphere.forEachParameters(
	// 	[&] (ECFParameter* p) { this->addParameter(p); }
	// );
}


bool EvolutionaryOptimizer::set(std::function<bool(void)> fnc)
{
	this->m_proxy.setNextConfiguration();
	this->m_set_function = fnc;
	
	int l_ret = fnc();
	int g_ret = 0;
	Communication::allReduce(&l_ret, &g_ret, 1, MPI_INT, MPI_MIN);
	
	if (!g_ret)
	{
		this->m_proxy.setConfigurationForbidden();
	}

	return static_cast<bool>(g_ret);
	// for (auto p = _parameters.begin(); p != _parameters.end(); ++p) {
	// 	std::cout << (*p)->name << ": " << (*p)->getValue() << " ";
	// }
	// std::cout << std::endl;
}

bool EvolutionaryOptimizer::run(std::function<bool(void)> fnc)
{
	int l_ret;
	
	double start = eslog::time();
	l_ret = fnc();
	double end = eslog::time();
	
	int g_ret = 0;
	Communication::allReduce(&l_ret, &g_ret, 1, MPI_INT, MPI_MIN);

	if (g_ret) { this->m_proxy.setConfigurationEvaluation(end - start); }
	else 
	{
		this->m_proxy.setConfigurationForbidden();
		while(!this->set(m_set_function));
	}

	return static_cast<bool>(g_ret);

	// this->m_proxy.setConfigurationEvaluation(sphere.evaluate());
}

