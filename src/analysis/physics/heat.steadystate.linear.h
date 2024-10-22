
#ifndef SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_
#define SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_

#include "analysis/physics/physics.h"
#include "analysis/assembler/heattransfer.h"
#include "analysis/pattern/pattern.h"
#include "analysis/math/vector_distributed.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct HeatTransferConfiguration;
struct HeatTransferLoadStepConfiguration;

class HeatSteadyStateLinear: public Physics {

public:
    HeatSteadyStateLinear(HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration);
    ~HeatSteadyStateLinear();

    bool analyze(step::Step &step);
    bool run(step::Step &step, Physics *prev);

    step::Time time;
    HeatTransferConfiguration &settings;
    HeatTransferLoadStepConfiguration &configuration;

    HeatTransfer assembler;

    Matrix_Base<double> *K;
    Vector_Distributed<Vector_Dense, double> *f, *x;
    Vector_Distributed<Vector_Sparse, double> *dirichlet;

    Pattern<double> *pattern;
    LinearSystemSolver<double> *solver;

protected:
    void storeSystem(step::Step &step);
    void storeSolution(step::Step &step);
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_HEAT_STEADYSTATE_LINEAR_H_ */
