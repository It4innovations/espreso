
#include "autoopt.h"
#include "ecfinfo.h"

#include <cstddef>

namespace espreso {
namespace autoopt {
namespace solver {

AutoOptimizer* autoopt = NULL;

void init(LoadStepSolverConfiguration &loadStep)
{
	if (autoopt != NULL) {
		delete autoopt;
	}
	if (info::ecf->auto_optimization.algorithm != AutoOptimizationConfiguration::ALGORITHM::NONE)
	{
		std::vector<ECFParameter*> opt_parameters;
		opt_parameters = {
			loadStep.feti.ecfdescription->getParameter(&loadStep.feti.preconditioner),
			loadStep.feti.ecfdescription->getParameter(&loadStep.feti.iterative_solver),
			loadStep.feti.ecfdescription->getParameter(&loadStep.feti.regularization),
			loadStep.feti.ecfdescription->getParameter(&loadStep.feti.redundant_lagrange),
			loadStep.feti.ecfdescription->getParameter(&loadStep.feti.B0_type),
			loadStep.feti.ecfdescription->getParameter(&loadStep.feti.scaling),
//			loadStep.feti.ecfdescription->getParameter(&loadStep.feti.method)
		};
		autoopt = new EvolutionaryOptimizer(info::ecf->auto_optimization, opt_parameters);
	} else {
		autoopt = new EmptyOptimizer();
	}
}

void update(std::function<bool(void)> fnc)
{
	while (!autoopt->set(fnc));
}

bool evaluate(std::function<bool(void)> fnc)
{
	return autoopt->run(fnc);
}

}
}
}
