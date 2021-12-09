
#ifndef SRC_AUTOOPT_OPTIMIZER_H_
#define SRC_AUTOOPT_OPTIMIZER_H_

#include "proxy.h"
#include "config/ecf/autoopt.h"

#include <vector>
#include <functional>

namespace espreso {

class AutoOptimizer {

public:
	virtual ~AutoOptimizer() {}

	virtual bool call(std::function<bool(void)> fnc)=0;

protected:
	AutoOptimizer() {}
};

class EmptyOptimizer : public AutoOptimizer {

public:
	EmptyOptimizer() {}

	bool call(std::function<bool(void)> fnc) override;
};

class EvolutionaryOptimizer : public AutoOptimizer {

public:
	EvolutionaryOptimizer(const AutoOptimizationConfiguration& configuration, std::vector<ECFParameter*>& parameters);

	bool call(std::function<bool(void)> fnc) override;

protected:
	OptimizationProxy m_proxy;
};

}



#endif /* SRC_AUTOOPT_OPTIMIZER_H_ */
