
#ifndef SRC_ASSEMBLER_INSTANCE_LINEAR_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_LINEAR_INSTANCE_H_

#include "../instance.h"

#include "../../../config/output.h"
#include "../../constraints/constraints.h"
#include "../../../output/vtk/vtk.h"
#include "../../../solver/generic/LinearSolver.h"

namespace espreso {

class ESPRESOSolver;

template <class TPhysics, class TConfiguration>
struct LinearInstance: public Instance
{
	LinearInstance(const TConfiguration &configuration, const OutputConfiguration &output, Mesh &mesh): Instance(mesh),
	_output(output),
	_configuration(configuration.espreso),
	_constrains(configuration.espreso, mesh),
	_physics(mesh, _constrains, configuration),
	_linearSolver(configuration.espreso, _physics, _constrains),
	_store(_output, mesh, "results")
	{
		_timeStatistics.totalTime.startWithBarrier();
	};

	virtual void init();
	virtual void solve(std::vector<std::vector<double> > &solution);
	virtual void finalize();

	virtual ~LinearInstance() {};

	virtual const OldPhysics& physics() const { return _physics; }
	virtual const Constraints& constraints() const { return _constrains; }
	virtual OldPhysics& physics() { return _physics; }
	virtual Constraints& constraints() { return _constrains; }

protected:
	const OutputConfiguration &_output;
	const ESPRESOSolver &_configuration;
	Constraints _constrains;
	TPhysics _physics;
	LinearSolver _linearSolver;
	store::VTK _store;
};

}

#include "instance.hpp"

#endif /* SRC_ASSEMBLER_INSTANCE_LINEAR_INSTANCE_H_ */
