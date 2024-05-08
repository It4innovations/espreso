
#ifndef SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_HARMONIC_REAL_LINEAR_H_
#define SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_HARMONIC_REAL_LINEAR_H_

#include "analysis/physics/physics.h"
#include "analysis/assembler/structuralmechanics.h"
#include "analysis/builder/builder.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsHarmonicRealLinear: public Physics {

public:
    StructuralMechanicsHarmonicRealLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration);
    ~StructuralMechanicsHarmonicRealLinear();

    void analyze(step::Step &step);
    void run(step::Step &step);

    step::Frequency frequency;
    StructuralMechanicsConfiguration &settings;
    StructuralMechanicsLoadStepConfiguration &configuration;

    StructuralMechanics assembler;

    Matrix_Base<double> *K, *M, *C;
    struct {
        Vector_Base<double> *f, *x, *dirichlet;
    } re, im;

    SparseMatrixBuilder<double> *builderAssembler, *builderSolver;
    LinearSystemSolver<double> *solver;

protected:
    void storeSystem(step::Step &step);
    void storeSolution(step::Step &step);
};

}

#endif /* SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_HARMONIC_REAL_LINEAR_H_ */
