
#ifndef SRC_PHYSICS_SYSTEM_SUPERLUSYSTEM_H_
#define SRC_PHYSICS_SYSTEM_SUPERLUSYSTEM_H_

#include "wrappers/superlu/w.superlu.systemsolver.h"
#include "distributedsystem.h"

namespace espreso {

struct SuperLUSolverData: public DistributedSolverData {

	SuperLUSolverData(SuperLUConfiguration &configuration)
	: DistributedSolverData(&solver),
	  solver(configuration, *this) {};

	SuperLUSystemSolver solver;
};

struct SuperLUSystem: public LinearSystem {

	std::vector<DistributedAssemblerData> assemblers;
	std::vector<SuperLUSolverData> solvers;

	virtual int nassemblers() { return assemblers.size(); }
	virtual int nsolvers() { return solvers.size(); }

	DistributedAssemblerData* assembler(int index = 0) { return assemblers.data() + index; }
	SuperLUSolverData* solver(int index = 0) { return solvers.data() + index; }

	SuperLUSystem(int assemblers, int solvers, SuperLUConfiguration &configuration);

protected:
	void _builderInit();
	void _builderReset();
	void _builderCreateSystem();
	void _builderUpdateSolution();
};

}

#endif /* SRC_PHYSICS_SYSTEM_SUPERLUSYSTEM_H_ */
