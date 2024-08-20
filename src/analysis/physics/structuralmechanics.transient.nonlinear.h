
#ifndef SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_TRANSIENT_NONLINEAR_H_
#define SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_TRANSIENT_NONLINEAR_H_

#include "analysis/physics/physics.h"
#include "analysis/assembler/structuralmechanics.h"
#include "analysis/builder/builder.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsTransientNonLinear: public Physics {

public:
    StructuralMechanicsTransientNonLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration);
    ~StructuralMechanicsTransientNonLinear();

    bool analyze(step::Step &step);
    bool run(step::Step &step);

    step::Time time;
    StructuralMechanicsConfiguration &settings;
    StructuralMechanicsLoadStepConfiguration &configuration;

    StructuralMechanics assembler;

    Matrix_Base<double> *K, *M, *C;
    Vector_Base<double> *f, *x, *dirichlet, *prev;
    Vector_Base<double> *U, *R, *dU, *V, *W, *X, *Y, *Z, *dTK, *dTM;

    SparseMatrixBuilder<double> *builder;
    LinearSystemSolver<double> *solver;

protected:
    bool checkDisplacement(step::Step &step);
    void storeSystem(step::Step &step);
    void storeSolution(step::Step &step);
};

}

#endif /* SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_TRANSIENT_NONLINEAR_H_ */
