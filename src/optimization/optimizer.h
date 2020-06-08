
#ifndef SRC_OPTIMIZATION_OPTIMIZER_H_
#define SRC_OPTIMIZATION_OPTIMIZER_H_

#include "../config/configuration.h"
#include "../config/ecf/solver/optimization/optimization.h"
#include "proxy.h"
#include "sphereproblem.h"

#include <vector>
#include <functional>

namespace espreso {

class Optimizer {

public:
	virtual ~Optimizer() {}

	virtual void set()=0;
	virtual void run(std::function<void(void)> fnc)=0;

protected:
	Optimizer() {}
};

class EmptyOptimizer : public Optimizer {

public:
	EmptyOptimizer() {}

	void set() override {}
	void run(std::function<void(void)> fnc) override;
};

class EvolutionaryOptimizer : public Optimizer {

public:
	EvolutionaryOptimizer(const OptimizationConfiguration& configuration,
		std::vector<ECFParameter*>& parameters);

	void set() override;
	void run(std::function<void(void)> fnc) override;

protected:
	OptimizationProxy proxy;
	SphereProblem sphere;
};

}



#endif /* SRC_OPTIMIZATION_OPTIMIZER_H_ */
