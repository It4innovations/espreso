
#ifndef SRC_PHYSICS_SYSTEM_PARDISOSYSTEM_H_
#define SRC_PHYSICS_SYSTEM_PARDISOSYSTEM_H_

#include "wrappers/pardiso/w.pardiso.systemsolver.h"
#include "distributedsystem.h"

namespace espreso {

struct PARDISOSolverData: public DistributedSolverData {

	PARDISOSolverData(PARDISOConfiguration &configuration)
	: DistributedSolverData(&solver),
	  solver(configuration, *this) {};

	PARDISOSystemSolver solver;
};

struct PARDISOSystem: public LinearSystem {

	std::vector<DistributedAssemblerData> assemblers;
	std::vector<PARDISOSolverData> solvers;

	virtual int nassemblers() { return assemblers.size(); }
	virtual int nsolvers() { return solvers.size(); }

	DistributedAssemblerData* assembler(int index = 0) { return assemblers.data() + index; }
	PARDISOSolverData* solver(int index = 0) { return solvers.data() + index; }

	PARDISOSystem(int assemblers, int solvers, PARDISOConfiguration &configuration);

protected:
	void _builderInit();
	void _builderReset();
	void _builderCreateSystem();
	void _builderUpdateSolution();
};

}


#endif /* SRC_PHYSICS_SYSTEM_PARDISOSYSTEM_H_ */
