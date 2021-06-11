
#ifndef SRC_ESINFO_AUTOOPT_H_
#define SRC_ESINFO_AUTOOPT_H_

#include "autoopt/optimizer.h"

#include <functional>

namespace espreso {

struct LoadStepSolverConfiguration;
struct DecompositionConfiguration;

namespace autoopt {
namespace solver {

	void init(LoadStepSolverConfiguration &loadStep, 
		DecompositionConfiguration &decomposition);

	void update(std::function<bool(void)> fnc);
	bool evaluate(std::function<bool(void)> fnc);

	extern AutoOptimizer *optimizer;
};

}
}




#endif /* SRC_ESINFO_AUTOOPT_H_ */
