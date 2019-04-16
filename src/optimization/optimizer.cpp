
#include "optimizer.h"

#include "../basis/logging/logging.h"

using namespace espreso;

Optimizer::Optimizer() 
: proxy(_parameters, OptimizationAlgorithm::SOMAT3A)
{
	// sphere.forEachParameters(
	// 	[&] (ECFParameter* p) { this->addParameter(p); }
	// );
}


void Optimizer::set()
{
	this->proxy.setNextConfiguration();
	// for (auto p = _parameters.begin(); p != _parameters.end(); ++p) {
	// 	std::cout << (*p)->name << ": " << (*p)->getValue() << " ";
	// }
	// std::cout << std::endl;
}

void Optimizer::run(std::function<void(void)> fnc)
{
	// double start = Measure::time();
	// fnc();
	// double end = Measure::time();

	// this->proxy.setConfigurationEvaluation(end - start);

	// this->proxy.setConfigurationEvaluation(sphere.evaluate());
}


