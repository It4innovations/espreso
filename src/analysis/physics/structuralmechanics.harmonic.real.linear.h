
#ifndef SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_HARMONIC_REAL_LINEAR_H_
#define SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_HARMONIC_REAL_LINEAR_H_

#include "analysis/physics/physics.h"
#include "analysis/assembler/structuralmechanics.h"
#include "analysis/pattern/pattern.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsHarmonicRealLinear: public Physics {

public:
    StructuralMechanicsHarmonicRealLinear(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration);
    ~StructuralMechanicsHarmonicRealLinear();

    bool analyze(step::Step &step);
    bool run(step::Step &step, Physics *prev);

    step::Frequency frequency;
    StructuralMechanicsConfiguration &settings;
    StructuralMechanicsLoadStepConfiguration &configuration;

    StructuralMechanics assembler;

    Matrix_Base<double> *K, *M, *C;
    struct {
        Vector_Distributed<Vector_Dense, double> *f, *x;
        Vector_Distributed<Vector_Sparse, double> *dirichlet;
    } re, im;

    Pattern<double> *patternAssembler, *patternSolver;
    LinearSystemSolver<double> *solver;

protected:
    void storeSystem(step::Step &step);
    void storeSolution(step::Step &step);
};

}

#endif /* SRC_ANALYSIS_PHYSICS_STRUCTURALMECHANICS_HARMONIC_REAL_LINEAR_H_ */
