
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
	HypreInstance(Mesh &mesh): Instance(mesh),
	feiPtr(MPI_COMM_WORLD),
	_constrains(mesh),
	_physics(mesh, _constrains),
	_store(mesh, "results", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO)
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
	LLNL_FEI_Impl feiPtr;
	Constraints _constrains;
	TPhysics _physics;
	output::VTK _store;

};

}

#include "../hypre/instance.hpp"

#else

namespace espreso {

template <class TPhysics>
struct HypreInstance: public Instance
{
public:
	HypreInstance(Mesh &mesh): Instance(mesh), _constrains(mesh), _physics(mesh, _constrains)
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
	Constraints _constrains;
	TPhysics _physics;
};

}

#endif // HAVE_HYPRE

#endif /* SRC_ASSEMBLER_INSTANCE_HYPRE_INSTANCE_H_ */
