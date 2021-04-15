
#include "optimizer.h"
#include "esinfo/eslog.h"

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
	
	bool ret = fnc();

	if (!ret)
	{
		this->m_proxy.setConfigurationForbidden();
	}

	return ret;
	// for (auto p = _parameters.begin(); p != _parameters.end(); ++p) {
	// 	std::cout << (*p)->name << ": " << (*p)->getValue() << " ";
	// }
	// std::cout << std::endl;
}

bool EvolutionaryOptimizer::run(std::function<bool(void)> fnc)
{
	bool ret;
	
	double start = eslog::time();
	ret = fnc();
	double end = eslog::time();

	if (ret) { this->m_proxy.setConfigurationEvaluation(end - start); }
	else 
	{
		this->m_proxy.setConfigurationForbidden();
		while(!this->set(m_set_function));
	}

	return ret;

	// this->m_proxy.setConfigurationEvaluation(sphere.evaluate());
}

