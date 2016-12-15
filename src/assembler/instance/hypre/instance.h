
#ifndef SRC_ASSEMBLER_INSTANCE_HYPRE_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_HYPRE_INSTANCE_H_

#include "../instance.h"
#include "../../../config/solverhypre.h"

#ifdef HAVE_HYPRE

#include "LLNL_FEI_Impl.h"

namespace espreso {

template <class TPhysics, class TConfiguration>
struct HypreInstance: public Instance
{
public:
	HypreInstance(const TConfiguration &configuration, const OutputConfiguration &output, Mesh &mesh): Instance(mesh),
	_output(output),
	_configuration(configuration.hypre),
	feiPtr(MPI_COMM_WORLD),
	_constrains(configuration.espreso, mesh),
	_physics(mesh, _constrains, configuration),
	_store(output, mesh, "results")
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
	const OutputConfiguration &_output;
	const HypreSolver &_configuration;
	LLNL_FEI_Impl feiPtr;
	Constraints _constrains;
	TPhysics _physics;
	store::VTK _store;

};

}

#include "../hypre/instance.hpp"

#else

namespace espreso {

template <class TPhysics, class TConfiguration>
struct HypreInstance: public Instance
{
public:
	HypreInstance(const TConfiguration &configuration, const OutputConfiguration &output, Mesh &mesh)
	: Instance(mesh), _configuration(configuration.hypre), _constrains(configuration.espreso, mesh), _physics(mesh, _constrains, configuration)
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
	Constraints _constrains;
	TPhysics _physics;
};

}

#endif // HAVE_HYPRE

#endif /* SRC_ASSEMBLER_INSTANCE_HYPRE_INSTANCE_H_ */
