
#ifndef SRC_PHYSICS_LOADSTEPSOLVER_TOPOLOGYOPTIMIZATION_H_
#define SRC_PHYSICS_LOADSTEPSOLVER_TOPOLOGYOPTIMIZATION_H_

#include "loadstepsolver.h"

namespace espreso {

struct TopologyOptimizationConfiguration;
class VectorDense;

class TopologyOptimization: public LoadStepSolver {

public:
	TopologyOptimization(LinearSystem *system, SubStepSolver *subStepSolver, TopologyOptimizationConfiguration &configuration);
	~TopologyOptimization();

	void init(LoadStepSolver *previous);
	void updateSystem();
	void updateStructuralMatrices();

protected:
	bool runNextSubstep() override;

	TopologyOptimizationConfiguration &_configuration;

	VectorDense *xPhys, *x, *DC, *C, *DV;
};

}

#endif /* SRC_PHYSICS_LOADSTEPSOLVER_TOPOLOGYOPTIMIZATION_H_ */
