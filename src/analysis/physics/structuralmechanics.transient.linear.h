
#ifndef SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_TRANSIENT_LINEAR_H_
#define SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_TRANSIENT_LINEAR_H_

#include "analysis/physics/physics.h"
#include "analysis/assembler/structuralmechanics.h"
#include "analysis/builder/builder.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsTransientLinear: public Physics {

public:
    StructuralMechanicsTransientLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration);
    ~StructuralMechanicsTransientLinear();

    bool analyze(step::Step &step);
    bool run(step::Step &step, Physics *prev);

    step::Time time;
    StructuralMechanicsConfiguration &settings;
    StructuralMechanicsLoadStepConfiguration &configuration;

    StructuralMechanics assembler;

    Matrix_Base<double> *K, *M;
    Vector_Base<double> *f, *x, *dirichlet, *checkpoint;
    Vector_Base<double> *U, *dU, *V, *W, *X, *Y, *Z, *dTK, *dTM;

    SparseMatrixBuilder<double> *builder;
    LinearSystemSolver<double> *solver;

protected:
    void storeSystem(step::Step &step);
    void storeSolution(step::Step &step);
};

}


#endif // SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_TRANSIENT_LINEAR_H_
