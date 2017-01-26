
#ifndef SRC_ASSEMBLER_INSTANCE_PRECOMPUTED_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_PRECOMPUTED_INSTANCE_H_

#include "../instance.h"

#include "../../../mesh/structures/mesh.h"
#include "../../../mesh/elements/element.h"
#include "../../../solver/generic/LinearSolver.h"

namespace espreso {

class ESPRESOSolver;

template <class TPhysics>
struct PrecomputedInstance: public Instance
{
	PrecomputedInstance(const ESPRESOSolver &configuration, APIMesh &mesh, SparseMatrix::MatrixType type, double* rhs, eslocal rhs_size)
	: Instance(mesh),
	  _configuration(configuration),
	  _constrains(configuration, mesh),
	  _physics(mesh, _constrains, configuration, type, rhs, rhs_size),
	  _linearSolver(configuration, _physics, _constrains)
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
	const ESPRESOSolver &_configuration;
	Constraints _constrains;
	TPhysics _physics;
	LinearSolver _linearSolver;
};

}

#include "../precomputed/instance.hpp"


#endif /* SRC_ASSEMBLER_INSTANCE_PRECOMPUTED_INSTANCE_H_ */
