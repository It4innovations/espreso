
#ifndef SRC_OPTIMIZATION_OPTIMIZER_H_
#define SRC_OPTIMIZATION_OPTIMIZER_H_

#include "../config/configuration.h"
#include "proxy.h"
#include "sphereproblem.h"

#include <vector>
#include <functional>

namespace espreso {

class Optimizer {

public:
	Optimizer();

	void addParameter(ECFParameter* parameter)
	{
		_parameters.push_back(parameter);
	}

	void set();
	void run(std::function<void(void)> fnc);

protected:
	std::vector<ECFParameter*> _parameters;
	OptimizationProxy proxy;
	SphereProblem sphere;
};

}



#endif /* SRC_OPTIMIZATION_OPTIMIZER_H_ */
