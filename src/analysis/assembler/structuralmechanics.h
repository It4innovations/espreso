
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_H_

#include "analysis/math/matrix_base.h"
#include "analysis/math/vector_base.h"
#include "analysis/math/vector_distributed.h"
#include "assembler.h"
#include "analysis/assembler/structuralmechanics/operators.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "config/holders/expression.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_info.h"

#include <cstddef>
#include <map>

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanicsLoadStepConfiguration;
struct SteadyState;

class StructuralMechanics: public Assembler
{
public:
    StructuralMechanics(StructuralMechanics *previous, StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration);

    bool analyze(const step::Step &step);
    void getInitialVelocity(Vector_Base<double> *x);

    void connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet);
    void evaluate(const step::Step &step, const step::Time &time, Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet);
    void updateSolution(Vector_Distributed<Vector_Dense, double> *x);
    void nextIteration(Vector_Distributed<Vector_Dense, double> *x);

    void connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Matrix_Base<double> *C, Vector_Base<double> *ref, Vector_Base<double> *imf, Vector_Base<double> *renf, Vector_Base<double> *imnf, Vector_Base<double> *reDirichlet, Vector_Base<double> *imDirichlet);
    void evaluate(const step::Step &step, const step::Frequency &freq, Matrix_Base<double> *K, Matrix_Base<double> *M, Matrix_Base<double> *C, Vector_Base<double> *ref, Vector_Base<double> *imf, Vector_Base<double> *renf, Vector_Base<double> *imnf, Vector_Base<double> *reDirichlet, Vector_Base<double> *imDirichlet);
    void updateSolution(Vector_Base<double> *rex, Vector_Base<double> *imx);

    StructuralMechanicsConfiguration &settings;
    StructuralMechanicsLoadStepConfiguration &configuration;

    struct Results {
        static NodeData *thickness;
        static NodeData *normal;
        static NodeData *initialVelocity;

        // steady state, transient
        static NodeData *displacement;
        static ElementData *principalStress, *componentStress, *vonMisesStress;
        static ElementData *isPlastized;

        // harmonic
        static NodeData *cosDisplacement, *sinDisplacement, *displacementAmplitude;
        static NodeData *phase, *velocity, *velocityAmplitude, *acceleration, *accelerationAmplitude;
        static NodeData *reactionForce;

        // FSI
        static NodeData *fluidForce;
    };

protected:
    step::Step step;
    step::Time time;
    step::Frequency frequency;
    void elements(SubKernel::Action action, size_t interval);
    void boundary(SubKernel::Action action, size_t region, size_t interval);
    void nodes(SubKernel::Action action, size_t region, size_t interval);
    void bem(SubKernel::Action action, size_t domain, double *BETI);

    std::vector<StructuralMechanicsElementOperators> elementKernels;
    std::vector<std::vector<StructuralMechanicsFaceOperators> > faceKernels;
    std::vector<std::vector<StructuralMechanicsNodeOperators> > nodeKernels;
    std::vector<double> faceMultiplicity;

    Vector_FETI<Vector_Dense, double> xBEM;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_H_ */
