
#ifndef SRC_ASSEMBLER_INSTANCE_HYPRE_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_HYPRE_INSTANCE_H_

#include "../instance.h"
#include "LLNL_FEI_Impl.h"

namespace espreso {

template <class TConstrains, class TPhysics>
struct HypreInstance: public Instance
{
public:
	HypreInstance(Mesh &mesh): Instance(mesh), feiPtr(MPI_COMM_WORLD), _constrains(mesh), _physics(mesh, _constrains)
	{
		_timeStatistics.totalTime.startWithBarrier();
	}

	virtual void init();
	virtual void solve(std::vector<std::vector<double> > &solution);
	virtual void finalize() {};

	virtual ~HypreInstance() {};

protected:
	LLNL_FEI_Impl feiPtr;
	TConstrains _constrains;
	TPhysics _physics;

};

}

#include "../hypre/instance.hpp"



#endif /* SRC_ASSEMBLER_INSTANCE_HYPRE_INSTANCE_H_ */
