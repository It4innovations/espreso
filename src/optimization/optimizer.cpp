
#include "optimizer.h"

#include "../basis/logging/logging.h"

using namespace espreso;

void Optimizer::set()
{
	for (auto p = _parameters.begin(); p != _parameters.end(); ++p) {
		std::cout << (*p)->name << ": " << (*p)->getValue() << "\n";
	}
}

void Optimizer::run(std::function<void(void)> fnc)
{
	double start = Measure::time();
	fnc();
	double end = Measure::time();
}


