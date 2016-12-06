
#ifndef SRC_ASSEMBLER_INSTANCE_SSNM_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_SSNM_INSTANCE_H_

#include "../../instance.h"
#include "esoutput.h"

namespace espreso {

template <class TPhysics>
struct SemiSmoothNewtonMethod: public Instance
{
	SemiSmoothNewtonMethod(Mesh &mesh): Instance(mesh),
	_constrains(mesh),
	_physics(mesh, _constrains),
	_linearSolver(_physics, _constrains),
	_store(mesh, "results", output->domain_shrink_ratio, output->cluster_shrink_ratio)
	{
		_timeStatistics.totalTime.startWithBarrier();
	};

	virtual void init();
	virtual void solve(std::vector<std::vector<double> > &solution);
	virtual void finalize();

	virtual ~SemiSmoothNewtonMethod() {};

	virtual const Physics& physics() const { return _physics; }
	virtual const Constraints& constraints() const { return _constrains; }

protected:
	Constraints _constrains;
	TPhysics _physics;
	LinearSolver _linearSolver;
	store::VTK _store;
};

}

#include "instance.hpp"

#endif /* SRC_ASSEMBLER_INSTANCE_SSNM_INSTANCE_H_ */
