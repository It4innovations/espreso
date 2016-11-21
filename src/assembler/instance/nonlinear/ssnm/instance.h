
#ifndef SRC_ASSEMBLER_INSTANCE_SSNM_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_SSNM_INSTANCE_H_

#include "../../instance.h"
#include "esoutput.h"
#include "../../../../config/output.h"

namespace espreso {

template <class TPhysics, class TConfiguration>
struct SemiSmoothNewtonMethod: public Instance
{
	SemiSmoothNewtonMethod(const TConfiguration &configuration, const OutputConfiguration &output, Mesh &mesh): Instance(mesh),
	_output(output),
	_configuration(configuration.espreso),
	_constraints(mesh),
	_physics(mesh, _constraints),
	_linearSolver(configuration.espreso, _physics, _constraints),
	_prevSolution(mesh.parts()),
	_store(_output, mesh, "results")
	{
		_timeStatistics.totalTime.startWithBarrier();
	};

	virtual void init();
	virtual void solve(std::vector<std::vector<double> > &solution);
	virtual void finalize();

	virtual ~SemiSmoothNewtonMethod() {};

	virtual const Physics& physics() const { return _physics; }
	virtual const Constraints& constraints() const { return _constraints; }

protected:
	const OutputConfiguration &_output;
	const ESPRESOSolver &_configuration;
	Constraints _constraints;
	TPhysics _physics;
	LinearSolver _linearSolver;
	std::vector<std::vector<double> > _prevSolution;
	store::VTK _store;
};

}

#include "instance.hpp"

#endif /* SRC_ASSEMBLER_INSTANCE_SSNM_INSTANCE_H_ */
