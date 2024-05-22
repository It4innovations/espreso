
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_H_

#include "analysis/math/matrix_base.h"
#include "analysis/math/vector_base.h"
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

    void analyze(const step::Step &step);

    void connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet);
    void evaluate(const step::Step &step, const step::Time &time, Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet);
    void updateSolution(const step::Step &step, Vector_Base<double> *x);
    void nextIteration(const step::Step &step, Vector_Base<double> *x);

    void connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Matrix_Base<double> *C, Vector_Base<double> *ref, Vector_Base<double> *imf, Vector_Base<double> *renf, Vector_Base<double> *imnf, Vector_Base<double> *reDirichlet, Vector_Base<double> *imDirichlet);
    void evaluate(const step::Step &step, const step::Frequency &freq, Matrix_Base<double> *K, Matrix_Base<double> *M, Matrix_Base<double> *C, Vector_Base<double> *ref, Vector_Base<double> *imf, Vector_Base<double> *renf, Vector_Base<double> *imnf, Vector_Base<double> *reDirichlet, Vector_Base<double> *imDirichlet);
    void updateSolution(const step::Step &step, Vector_Base<double> *rex, Vector_Base<double> *imx);

    StructuralMechanicsConfiguration &settings;
    StructuralMechanicsLoadStepConfiguration &configuration;

    struct Results {
        static NodeData *thickness;
        static NodeData *normal;

        // steady state, transient
        static NodeData *displacement;
        static ElementData *principalStress, *componentStress, *vonMisesStress;
        static ElementData *isPlastized;

        // harmonic
        static NodeData *cosDisplacement, *sinDisplacement, *displacementAmplitude;
        static NodeData *phase, *velocity, *velocityAmplitude, *acceleration, *accelerationAmplitude;
    };

    struct {
        bool K, M, C, f, nf, dirichlet;
    } constant;

protected:
    void run(const step::Step &step, SubKernel::Action action, size_t interval);
    void run(const step::Step &step, SubKernel::Action action, size_t region, size_t interval);
    void runBEM(const step::Step &step, SubKernel::Action action, size_t domain, double *BETI);

    std::vector<StructuralMechanicsElementOperators> subkernels;
    std::vector<std::vector<StructuralMechanicsBoundaryOperators> > boundary;
    std::vector<double> faceMultiplicity;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_H_ */
