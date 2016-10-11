
#ifndef SRC_ASSEMBLER_INSTANCE_PRECOMPUTED_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_PRECOMPUTED_INSTANCE_H_

#include "../instance.h"

namespace espreso {

template <class TPhysics>
struct PrecomputedInstance: public Instance
{
	PrecomputedInstance(APIMesh &mesh, SparseMatrix::MatrixType type, double* rhs, eslocal rhs_size)
	: Instance(mesh), _constrains(mesh), _physics(mesh, _constrains, type, rhs, rhs_size), _linearSolver(_physics, _constrains)
	{
		_timeStatistics.totalTime.startWithBarrier();
	};

	virtual void init();
	virtual void solve(std::vector<std::vector<double> > &solution);
	virtual void finalize();

	virtual ~PrecomputedInstance() {};

	virtual const Physics& physics() const { return _physics; }
	virtual const Constraints& constraints() const { return _constrains; }

protected:
	Constraints _constrains;
	TPhysics _physics;
	LinearSolver _linearSolver;
};

}

#include "../precomputed/instance.hpp"


#endif /* SRC_ASSEMBLER_INSTANCE_PRECOMPUTED_INSTANCE_H_ */
