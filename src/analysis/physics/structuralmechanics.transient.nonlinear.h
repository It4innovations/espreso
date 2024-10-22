
#ifndef SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_TRANSIENT_NONLINEAR_H_
#define SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_TRANSIENT_NONLINEAR_H_

#include "analysis/physics/physics.h"
#include "analysis/assembler/structuralmechanics.h"
#include "analysis/pattern/pattern.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsTransientNonLinear: public Physics {

public:
    StructuralMechanicsTransientNonLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration);
    ~StructuralMechanicsTransientNonLinear();

    bool analyze(step::Step &step);
    bool run(step::Step &step, Physics *prev);

    step::Time time;
    StructuralMechanicsConfiguration &settings;
    StructuralMechanicsLoadStepConfiguration &configuration;

    StructuralMechanics assembler;

    Matrix_Base<double> *K, *M, *C;
    Vector_Distributed<Vector_Dense, double> *f, *f_old;
    Vector_Distributed<Vector_Dense, double> *R, *R_old, *dU, *U, *V, *A, *U_old, *V_old, *A_old, *X;
    Vector_Distributed<Vector_Sparse, double> *dirichlet;

    Pattern<double> *pattern;
    LinearSystemSolver<double> *solver;

protected:
    bool checkDisplacement(step::Step &step, double f_norm);
    void storeSystem(step::Step &step);
    void storeSolution(step::Step &step);
};

}

#endif /* SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_TRANSIENT_NONLINEAR_H_ */
