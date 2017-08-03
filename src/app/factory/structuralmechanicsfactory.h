
#ifndef SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_
#define SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_

#include "factory.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanics2DConfiguration;
struct StructuralMechanics3DConfiguration;;

class StructuralMechanicsFactory: public FactoryLoader {

public:
	StructuralMechanicsFactory(const StructuralMechanics2DConfiguration &configuration, Mesh *mesh);
	StructuralMechanicsFactory(const StructuralMechanics3DConfiguration &configuration, Mesh *mesh);

	size_t loadSteps() const;
	LoadStepSolver* getLoadStepSolver(size_t step, Mesh *mesh, output::Store *store);

protected:
	const StructuralMechanicsConfiguration &_configuration;
	bool _bem;

};

}


#endif /* SRC_APP_FACTORY_STRUCTURALMECHANICSFACTORY_H_ */
