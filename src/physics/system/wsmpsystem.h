
#ifndef SRC_PHYSICS_SYSTEM_WSMPSYSTEM_H_
#define SRC_PHYSICS_SYSTEM_WSMPSYSTEM_H_

#include "wrappers/wsmp/w.wsmp.systemsolver.h"
#include "distributedsystem.h"

namespace espreso {

struct WSMPSolverData: public DistributedSolverData {

	WSMPSolverData(WSMPConfiguration &configuration)
	: DistributedSolverData(&solver),
	  solver(configuration, *this) {};

	WSMPSystemSolver solver;
};

struct WSMPSystem: public LinearSystem {

	std::vector<DistributedAssemblerData> assemblers;
	std::vector<WSMPSolverData> solvers;

	virtual int nassemblers() { return assemblers.size(); }
	virtual int nsolvers() { return solvers.size(); }

	DistributedAssemblerData* assembler(int index = 0) { return assemblers.data() + index; }
	WSMPSolverData* solver(int index = 0) { return solvers.data() + index; }

	WSMPSystem(int assemblers, int solvers, WSMPConfiguration &configuration);

protected:
	void _builderInit();
	void _builderReset();
	void _builderCreateSystem();
	void _builderUpdateSolution();
};

}

#endif /* SRC_PHYSICS_SYSTEM_WSMPSYSTEM_H_ */
