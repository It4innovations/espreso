
#ifndef SRC_ASSEMBLER_INSTANCE_PRECOMPUTED_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_PRECOMPUTED_INSTANCE_H_

#include "../instance.h"

namespace espreso {

template <class TConstrains, class TPhysics>
struct PrecomputedInstance: public Instance
{
	PrecomputedInstance(APIMesh &mesh, double* rhs, eslocal rhs_size, eslocal dirichlet_size, eslocal* dirichlet_indices, double* dirichlet_values, eslocal indexBase)
	: Instance(mesh), _constrains(mesh /*, dirichlet_size, dirichlet_indices, dirichlet_values, indexBase*/), _physics(mesh, _constrains, rhs, rhs_size), _linearSolver(_physics, _constrains)
	{
		_timeStatistics.totalTime.startWithBarrier();
	};

	virtual void init();
	virtual void solve(std::vector<std::vector<double> > &solution);
	virtual void finalize();

	virtual ~PrecomputedInstance() {};

protected:
	TConstrains _constrains;
	TPhysics _physics;
	LinearSolver _linearSolver;
};

}

#include "../precomputed/instance.hpp"


#endif /* SRC_ASSEMBLER_INSTANCE_PRECOMPUTED_INSTANCE_H_ */
