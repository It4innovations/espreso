
#ifndef SRC_PHYSICS_SYSTEM_BUILDER_TIMEBUILDER_H_
#define SRC_PHYSICS_SYSTEM_BUILDER_TIMEBUILDER_H_

#include "builder.h"

namespace espreso {

struct TimeBuilder: public Builder {
	void init(WSMPSystem &system);
	void init(SuperLUSystem &system);
	void init(PARDISOSystem &system);
	void init(MKLPDSSSystem &system);
	void init(HYPRESystem &system);
	void init(FETISystem &system);

	void buildSystem(WSMPSystem &system);
	void buildSystem(SuperLUSystem &system);
	void buildSystem(PARDISOSystem &system);
	void buildSystem(MKLPDSSSystem &system);
	void buildSystem(HYPRESystem &system);
	void buildSystem(FETISystem &system);

	void updateSolution(WSMPSystem &system);
	void updateSolution(SuperLUSystem &system);
	void updateSolution(PARDISOSystem &system);
	void updateSolution(MKLPDSSSystem &system);
	void updateSolution(HYPRESystem &system);
	void updateSolution(FETISystem &system);

protected:
	void init(AssemblerData &assembler, SolverData &solver);
	void buildSystem(AssemblerData &assembler, SolverData &solver);
	void updateSolution(AssemblerData &assembler, SolverData &solver);
};

}

#endif /* SRC_PHYSICS_SYSTEM_BUILDER_TIMEBUILDER_H_ */
