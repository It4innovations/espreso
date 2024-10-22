
#ifndef SRC_ANALYSIS_ANALYSIS_PLASTICITY_STEADYSTATE_H_
#define SRC_ANALYSIS_ANALYSIS_PLASTICITY_STEADYSTATE_H_

#include "analysis/physics/physics.h"
#include "analysis/assembler/structuralmechanics.h"
#include "analysis/pattern/pattern.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsSteadyStateNonLinear: public Physics {

public:
    StructuralMechanicsSteadyStateNonLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration);
    ~StructuralMechanicsSteadyStateNonLinear();

    bool analyze(step::Step &step);
    bool run(step::Step &step, Physics *prev);

    step::Time time;
    StructuralMechanicsConfiguration &settings;
    StructuralMechanicsLoadStepConfiguration &configuration;

    StructuralMechanics assembler;

    Matrix_Base<double> *K;
    Vector_Distributed<Vector_Dense, double> *U, *R, *f, *x;
    Vector_Distributed<Vector_Sparse, double> *dirichlet;

    Pattern<double> *pattern;
    LinearSystemSolver<double> *solver;

protected:
    bool checkDisplacement(step::Step &step);

    void storeSystem(step::Step &step);
    void storeSolution(step::Step &step);
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_PLASTICITY_STEADYSTATE_H_ */
