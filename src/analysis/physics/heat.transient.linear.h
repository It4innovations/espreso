
#ifndef SRC_ANALYSIS_PHYSICS_HEAT_TRANSIENT_LINEAR_H_
#define SRC_ANALYSIS_PHYSICS_HEAT_TRANSIENT_LINEAR_H_

#include "analysis/physics/physics.h"
#include "analysis/assembler/heattransfer.h"
#include "analysis/pattern/pattern.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct HeatTransferConfiguration;
struct HeatTransferLoadStepConfiguration;

class HeatTransientLinear: public Physics {

public:
    HeatTransientLinear(HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration);
    ~HeatTransientLinear();

    bool analyze(step::Step &step);
    bool run(step::Step &step, Physics *prev);

    step::Time time;
    HeatTransferConfiguration &settings;
    HeatTransferLoadStepConfiguration &configuration;

    HeatTransfer assembler;

    Matrix_Base<double> *K, *M;
    Vector_Distributed<Vector_Dense, double> *f, *x;
    Vector_Distributed<Vector_Dense, double> *U, *dU, *V, *X, *Y, *dTK, *dTM;
    Vector_Distributed<Vector_Sparse, double> *dirichlet;

    Pattern<double> *pattern;
    LinearSystemSolver<double> *solver;

protected:
    void storeSystem(step::Step &step);
    void storeSolution(step::Step &step);
};

}



#endif /* SRC_ANALYSIS_PHYSICS_HEAT_TRANSIENT_LINEAR_H_ */
