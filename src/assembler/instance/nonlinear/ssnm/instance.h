
#ifndef SRC_ASSEMBLER_INSTANCE_SSNM_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_SSNM_INSTANCE_H_

#include "../../instance.h"
#include "esoutput.h"

namespace espreso {

template <class TPhysics>
struct SemiSmoothNewtonMethod: public Instance
{
	SemiSmoothNewtonMethod(const OutputConfiguration &output, Mesh &mesh): Instance(mesh),
	_output(output),
	_constrains(mesh),
	_physics(mesh, _constrains),
	_linearSolver(_physics, _constrains),
	_store(_output, mesh, "results")
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
	const OutputConfiguration &_output;
	Constraints _constrains;
	TPhysics _physics;
	LinearSolver _linearSolver;
	store::VTK _store;
};

}

#include "instance.hpp"

#endif /* SRC_ASSEMBLER_INSTANCE_SSNM_INSTANCE_H_ */
