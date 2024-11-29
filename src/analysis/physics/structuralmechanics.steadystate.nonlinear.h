
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

    Matrix_Base<double> *K, *postM;
    Vector_Distributed<Vector_Dense, double> *U, *R, *f;
    Vector_Distributed<Vector_Sparse, double> *dirichlet;
    Vector_Distributed<Matrix_Dense, double> *postB, *postX;

    Pattern<double> *pattern, *postPattern;
    LinearSystemSolver<double> *solver, *postSolver;

protected:
    bool checkDisplacement(step::Step &step, double U_norm);
    bool checkStress(step::Step &step, double f_norm);

    void storeSystem(step::Step &step);
    void storeSolution(step::Step &step);
    void storePostSystem(step::Step &step);
    void storePostSolution(step::Step &step);
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_PLASTICITY_STEADYSTATE_H_ */
