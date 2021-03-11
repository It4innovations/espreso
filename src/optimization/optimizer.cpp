
#include "optimizer.h"

#include "../basis/logging/logging.h"

using namespace espreso;

bool EmptyOptimizer::set(std::function<bool(void)> fnc)
{
	return fnc();
}

bool EmptyOptimizer::run(std::function<bool(void)> fnc)
{
	return fnc();
}

EvolutionaryOptimizer::EvolutionaryOptimizer(const OptimizationConfiguration& configuration,
	std::vector<ECFParameter*>& parameters) 
: proxy(parameters, configuration)
{
	// sphere.forEachParameters(
	// 	[&] (ECFParameter* p) { this->addParameter(p); }
	// );
}


bool EvolutionaryOptimizer::set(std::function<bool(void)> fnc)
{
	this->proxy.setNextConfiguration();

	return fnc();
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

	if (ret) { this->proxy.setConfigurationEvaluation(end - start); }
	else { this->proxy.setConfigurationForbidden(); }

	return ret;

	// this->proxy.setConfigurationEvaluation(sphere.evaluate());
}

