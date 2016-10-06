
#ifndef SRC_ASSEMBLER_INSTANCE_LINEAR_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_LINEAR_INSTANCE_H_

#include "../instance.h"
#include "esoutput.h"

namespace espreso {

template <class TConstrains, class TPhysics>
struct LinearInstance: public Instance
{
public:
	LinearInstance(Mesh &mesh): Instance(mesh),
	_constrains(mesh),
	_physics(mesh, _constrains),
	_linearSolver(_physics, _constrains),
	_store(mesh, "results", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO)
	{
		_timeStatistics.totalTime.startWithBarrier();
	};

	virtual void init();
	virtual void solve(std::vector<std::vector<double> > &solution);
	virtual void finalize();

	virtual ~LinearInstance() {};

	virtual const Physics& physics() const { return _physics; }
	virtual const Constraints& constraints() const { return _constrains; }

protected:
	TConstrains _constrains;
	TPhysics _physics;
	LinearSolver _linearSolver;
	output::VTK _store;
};

}

#include "instance.hpp"

#endif /* SRC_ASSEMBLER_INSTANCE_LINEAR_INSTANCE_H_ */
