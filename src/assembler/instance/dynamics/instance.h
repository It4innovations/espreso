
#ifndef SRC_ASSEMBLER_INSTANCE_DYNAMICS_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_DYNAMICS_INSTANCE_H_

#include "../instance.h"

#include "../../../configuration/output.h"
#include "../../../old_output/vtk/vtk.h"
#include "../../constraints/constraints.h"
#include "../../../solver/generic/LinearSolver.h"

namespace espreso {

template <class TPhysics>
struct DynamicsInstance: public OldInstance
{
public:
	DynamicsInstance(const OutputConfiguration &output, Mesh &mesh): OldInstance(mesh),
	_output(output),
	_constrains(mesh),
	_physics(mesh, _constrains),
	_linearSolver(_physics, _constrains),
	_store(output, mesh, "results"),
	 _time(0)
	{
		_timeStatistics.totalTime.startWithBarrier();
	};

	virtual void init();

	virtual void solve(std::vector<std::vector<double> > &solution);

	virtual void finalize();

	virtual ~DynamicsInstance() {};

	virtual const OldPhysics& physics() const { return _physics; }
	virtual const Constraints& constraints() const { return _constrains; }
	virtual OldPhysics& physics() { return _physics; }
	virtual Constraints& constraints() { return _constrains; }

protected:
	const OutputConfiguration &_output;
	Constraints _constrains;
	TPhysics _physics;
	LinearSolver _linearSolver;
	store::VTK _store;

private:
	size_t _time;

	std::vector<std::vector<double> > _u, _v, _w;    // old vectors
	std::vector<std::vector<double> > _vn, _wn; // new vectors
	std::vector<std::vector<double> > _b, _tmp;
};

}

#include "instance.hpp"



#endif /* SRC_ASSEMBLER_INSTANCE_DYNAMICS_INSTANCE_H_ */
