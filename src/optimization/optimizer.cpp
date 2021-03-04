
#include "optimizer.h"

#include "../basis/logging/logging.h"

using namespace espreso;

void EmptyOptimizer::run(std::function<void(void)> fnc)
{
	fnc();
}

EvolutionaryOptimizer::EvolutionaryOptimizer(const OptimizationConfiguration& configuration,
	std::vector<ECFParameter*>& parameters) 
: proxy(parameters, configuration)
{
	// sphere.forEachParameters(
	// 	[&] (ECFParameter* p) { this->addParameter(p); }
	// );
}


void EvolutionaryOptimizer::set()
{
	this->proxy.setNextConfiguration();
	// for (auto p = _parameters.begin(); p != _parameters.end(); ++p) {
	// 	std::cout << (*p)->name << ": " << (*p)->getValue() << " ";
	// }
	// std::cout << std::endl;
}

void EvolutionaryOptimizer::run(std::function<void(void)> fnc)
{
	double start = Measure::time();
	fnc();
	double end = Measure::time();

	this->proxy.setConfigurationEvaluation(end - start);

	// this->proxy.setConfigurationEvaluation(sphere.evaluate());
}

