
#ifndef SRC_PHYSICS_LINEARSYSTEM_SYSTEM_HYPRESYSTEM_H_
#define SRC_PHYSICS_LINEARSYSTEM_SYSTEM_HYPRESYSTEM_H_

#include "math/vector.dense.distributed.h"
#include "wrappers/hypre/w.hypre.systemsolver.h"
#include "distributedsystem.h"

namespace espreso {

struct HYPREAssemblerData: public DistributedAssemblerData {
	int numFnc;
	VectorsDenseDistributed N;

	HYPREAssemblerData(): numFnc(1) {}

	void print(const Builder *builder, const char* prefix, const char* suffix);
};

struct HYPRESolverData: public DistributedSolverData {

	HYPRESolverData(HYPREConfiguration &configuration)
	: DistributedSolverData(&solver),
	  numFnc(1),
	  solver(configuration, *this) {};

	int numFnc;
	VectorsDenseDistributed N;
	HYPRESystemSolver solver;

	void printData(const Builder *builder, const char* prefix);
};

struct HYPRESystem: public LinearSystem {

	std::vector<HYPREAssemblerData> assemblers;
	std::vector<HYPRESolverData> solvers;

	virtual int nassemblers() { return assemblers.size(); }
	virtual int nsolvers() { return solvers.size(); }

	DistributedAssemblerData* assembler(int index = 0) { return assemblers.data() + index; }
	HYPRESolverData* solver(int index = 0) { return solvers.data() + index; }

	HYPRESystem(int assemblers, int solvers, HYPREConfiguration &configuration);

protected:
	void _builderInit();
	void _builderReset();
	void _builderCreateSystem();
	void _builderUpdateSolution();
};

}



#endif /* SRC_PHYSICS_LINEARSYSTEM_SYSTEM_HYPRESYSTEM_H_ */
