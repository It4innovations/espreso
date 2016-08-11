
#ifndef SRC_ASSEMBLER_INSTANCE_DYNAMICS_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_DYNAMICS_INSTANCE_H_

#include "../instance.h"

namespace espreso {

template <class TConstrains, class TPhysics>
struct DynamicsInstance: public Instance
{
public:
	DynamicsInstance(Mesh &mesh): Instance(mesh), _constrains(mesh), _physics(mesh, _constrains), _linearSolver(_physics, _constrains), _time(0)
	{
		_timeStatistics.totalTime.startWithBarrier();
	};

	virtual void init();

	virtual void pre_solve_update(std::vector<std::vector<double> > &solution);
	virtual void solve(std::vector<std::vector<double> > &solution);
	virtual void post_solve_update(std::vector<std::vector<double> > &solution);

	virtual void finalize();


	virtual ~DynamicsInstance() {};

protected:
	TConstrains _constrains;
	TPhysics _physics;
	LinearSolver _linearSolver;

private:
	size_t _time;

	std::vector<std::vector<double> > _u, _v, _w;    // old vectors
	std::vector<std::vector<double> > _vn, _wn; // new vectors
	std::vector<std::vector<double> > _b, _tmp;
};

}

#include "instance.hpp"



#endif /* SRC_ASSEMBLER_INSTANCE_DYNAMICS_INSTANCE_H_ */
