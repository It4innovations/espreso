
#ifndef SRC_ASSEMBLER_INSTANCE_LINEAR_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_LINEAR_INSTANCE_H_

#include "../instance.h"

namespace espreso {

template <class TConstrains, class TPhysics>
struct LinearInstance: public Instance
{
public:
	LinearInstance(const Mesh &mesh): Instance(mesh), _physics(mesh), _constrains(mesh, _physics.DOFs), _linearSolver(_physics, _constrains)
	{
		_timeStatistics.totalTime.startWithBarrier();
	};

	virtual void init();
	virtual void solve(std::vector<std::vector<double> > &solution);
	virtual void finalize();

	virtual ~LinearInstance() {};

protected:
	TPhysics _physics;
	TConstrains _constrains;
	LinearSolver _linearSolver;
};

}

#include "instance.hpp"

#endif /* SRC_ASSEMBLER_INSTANCE_LINEAR_INSTANCE_H_ */
