
#ifndef SRC_PHYSICS_LINEARSYSTEM_SYSTEM_MKLPDSSSYSTEM_H_
#define SRC_PHYSICS_LINEARSYSTEM_SYSTEM_MKLPDSSSYSTEM_H_

#include "wrappers/mklpdss/w.mkl.pdss.systemsolver.h"
#include "distributedsystem.h"

namespace espreso {

struct MKLPDSSSolverData: public DistributedSolverData {

    MKLPDSSSolverData(MKLPDSSConfiguration &configuration)
    : DistributedSolverData(&solver),
      solver(configuration, *this) {};

    MKLPDSSSystemSolver solver;
};

struct MKLPDSSSystem: public LinearSystem {

    std::vector<DistributedAssemblerData> assemblers;
    std::vector<MKLPDSSSolverData> solvers;

    virtual int nassemblers() { return assemblers.size(); }
    virtual int nsolvers() { return solvers.size(); }

    DistributedAssemblerData* assembler(int index = 0) { return assemblers.data() + index; }
    MKLPDSSSolverData* solver(int index = 0) { return solvers.data() + index; }

    MKLPDSSSystem(int assemblers, int solvers, MKLPDSSConfiguration &configuration);

protected:
    void _builderInit();
    void _builderReset();
    void _builderCreateSystem();
    void _builderUpdateSolution();
};

}

#endif /* SRC_PHYSICS_LINEARSYSTEM_SYSTEM_MKLPDSSSYSTEM_H_ */
