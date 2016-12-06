
#ifndef SRC_ASSEMBLER_INSTANCE_HYPRE_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_HYPRE_INSTANCE_H_

#include "../instance.h"

#ifdef HAVE_HYPRE

#include "LLNL_FEI_Impl.h"

namespace espreso {

template <class TPhysics>
struct HypreInstance: public Instance
{
public:
	HypreInstance(const HypreSolver &configuration, Mesh &mesh): Instance(mesh),
	_configuration(configuration),
	feiPtr(MPI_COMM_WORLD),
	_constrains(_dummyESPRESOSolver, mesh),
	_physics(mesh, _constrains, _dummyESPRESOSolver),
	_store(mesh, "results", output->domain_shrink_ratio, output->cluster_shrink_ratio)
	{
		_timeStatistics.totalTime.startWithBarrier();
	}

	virtual void init();
	virtual void solve(std::vector<std::vector<double> > &solution);
	virtual void finalize() {};

	virtual ~HypreInstance() {};

	virtual const Physics& physics() const { return _physics; }
	virtual const Constraints& constraints() const { return _constrains; }

protected:
	const HypreSolver &_configuration;
	const ESPRESOSolver _dummyESPRESOSolver;
	LLNL_FEI_Impl feiPtr;
	Constraints _constrains;
	TPhysics _physics;
	store::VTK _store;

};

}

#include "../hypre/instance.hpp"

#else

namespace espreso {

template <class TPhysics>
struct HypreInstance: public Instance
{
public:
	HypreInstance(const HypreSolver &configuration, Mesh &mesh): Instance(mesh), _configuration(configuration), _constrains(_dummyESPRESOSolver, mesh), _physics(mesh, _constrains, _dummyESPRESOSolver)
	{
		ESINFO(GLOBAL_ERROR) << "HYPRE is not linked! Specify HYPRE::INCLUDE and HYPRE::LIBPATH";
	}

	virtual void init() {};
	virtual void solve(std::vector<std::vector<double> > &solution) {};
	virtual void finalize() {};

	virtual ~HypreInstance() {};

	virtual const Physics& physics() const { return _physics; }
	virtual const Constraints& constraints() const { return _constrains; }

protected:
	const HypreSolver &_configuration;
	const ESPRESOSolver _dummyESPRESOSolver;
	Constraints _constrains;
	TPhysics _physics;
};

}

#endif // HAVE_HYPRE

#endif /* SRC_ASSEMBLER_INSTANCE_HYPRE_INSTANCE_H_ */
