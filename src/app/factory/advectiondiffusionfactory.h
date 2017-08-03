
#ifndef SRC_APP_FACTORY_ADVECTIONDIFFUSIONFACTORY_H_
#define SRC_APP_FACTORY_ADVECTIONDIFFUSIONFACTORY_H_

#include "factory.h"

namespace espreso {

struct AdvectionDiffusionConfiguration;
struct AdvectionDiffusion2DConfiguration;
struct AdvectionDiffusion3DConfiguration;

class AdvectionDiffusionFactory: public FactoryLoader {

public:
	AdvectionDiffusionFactory(const AdvectionDiffusion2DConfiguration &configuration, Mesh *mesh);
	AdvectionDiffusionFactory(const AdvectionDiffusion3DConfiguration &configuration, Mesh *mesh);

	size_t loadSteps() const;
	LoadStepSolver* getLoadStepSolver(size_t step, Mesh *mesh, output::Store *store);

protected:
	const AdvectionDiffusionConfiguration &_configuration;
	bool _bem;

};

}



#endif /* SRC_APP_FACTORY_ADVECTIONDIFFUSIONFACTORY_H_ */
