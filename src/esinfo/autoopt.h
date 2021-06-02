
#ifndef SRC_ESINFO_AUTOOPT_H_
#define SRC_ESINFO_AUTOOPT_H_

#include "autoopt/optimizer.h"

#include <functional>

namespace espreso {

struct LoadStepSolverConfiguration;

namespace autoopt {
namespace solver {

	void init(LoadStepSolverConfiguration &loadStep);

	void update(std::function<bool(void)> fnc);
	void evaluate(std::function<bool(void)> fnc);

	extern AutoOptimizer *optimizer;
};

}
}




#endif /* SRC_ESINFO_AUTOOPT_H_ */
