
#ifndef SRC_ASSEMBLER_INSTANCE_HYPRE_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_HYPRE_INSTANCE_H_

#include "../instance.h"
#include "../../../configuration/solver/hypre.h"
#include "../../../output/resultstore/vtkxmlascii.h"

#ifdef HAVE_HYPRE

#include "LLNL_FEI_Impl.h"

namespace espreso {

template <class TPhysics, class TConfiguration>
struct HypreInstance: public OldInstance
{
public:
	HypreInstance(const TConfiguration &configuration, const OutputConfiguration &output, Mesh &mesh): OldInstance(mesh),
	_output(output),
	_configuration(configuration.physics_solver.load_steps_settings.at(1)->hypre),
	feiPtr(MPI_COMM_WORLD),
	_constrains(configuration.physics_solver.load_steps_settings.at(1)->espreso, mesh),
	_physics(mesh, _constrains, configuration),
	_store(output, &mesh, "results")
	{
		_timeStatistics.totalTime.startWithBarrier();
	}

	virtual void init();
	virtual void solve(std::vector<std::vector<double> > &solution);
	virtual void finalize() {};

	virtual ~HypreInstance() {};

	virtual const OldPhysics& physics() const { return _physics; }
	virtual const Constraints& constraints() const { return _constrains; }
	virtual OldPhysics& physics() { return _physics; }
	virtual Constraints& constraints() { return _constrains; }

protected:
	const OutputConfiguration &_output;
	const HypreSolver &_configuration;
	LLNL_FEI_Impl feiPtr;
	Constraints _constrains;
	TPhysics _physics;
	output::VTKXMLASCII _store;

};

}

#include "../hypre/instance.hpp"

#else

namespace espreso {

template <class TPhysics, class TConfiguration>
struct HypreInstance: public OldInstance
{
public:
	HypreInstance(const TConfiguration &configuration, const OutputConfiguration &output, Mesh &mesh)
	: OldInstance(mesh), _configuration(configuration.physics_solver.load_steps_settings.at(1)->hypre), _constrains(configuration.physics_solver.load_steps_settings.at(1)->espreso, mesh), _physics(mesh, _constrains, configuration)
	{
		ESINFO(GLOBAL_ERROR) << "HYPRE is not linked! Specify HYPRE::INCLUDE and HYPRE::LIBPATH";
	}

	virtual void init() {};
	virtual void solve(std::vector<std::vector<double> > &solution) {};
	virtual void finalize() {};

	virtual ~HypreInstance() {};

	virtual const OldPhysics& physics() const { return _physics; }
	virtual const Constraints& constraints() const { return _constrains; }
	virtual OldPhysics& physics() { return _physics; }
	virtual Constraints& constraints() { return _constrains; }

protected:
	const HypreSolver &_configuration;
	Constraints _constrains;
	TPhysics _physics;
};

}

#endif // HAVE_HYPRE

#endif /* SRC_ASSEMBLER_INSTANCE_HYPRE_INSTANCE_H_ */
