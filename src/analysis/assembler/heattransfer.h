
#ifndef SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_H_
#define SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_H_

#include "analysis/math/matrix_base.h"
#include "analysis/math/vector_base.h"
#include "assembler.h"
#include "analysis/assembler/heattransfer/operators.h"
#include "config/ecf/physics/heattransfer.h"
#include "config/holders/expression.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_info.h"

#include <cstddef>
#include <map>

#include <type_traits>

namespace espreso {

struct HeatTransferConfiguration;
struct HeatTransferLoadStepConfiguration;
struct SteadyState;

class HeatTransfer: public Assembler
{
public:
    HeatTransfer(HeatTransfer *previous, HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration);

    void analyze();
    void getInitialTemperature(Vector_Base<double> *x);

    void connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet);
    void evaluate(step::Time &time, Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet);
    void updateSolution(Vector_Base<double> *x);

    HeatTransferConfiguration &settings;
    HeatTransferLoadStepConfiguration &configuration;

    struct Results {
        static NodeData *temperature, *initialTemperature, *thickness;
        static ElementData *translationMotion, *gradient, *flux;
    };

    struct {
        bool K, M, f, nf, dirichlet;
    } constant;

protected:
    void run(SubKernel::Action action, size_t interval);
    void run(SubKernel::Action action, size_t region, size_t interval);
    void runBEM(SubKernel::Action action, size_t domain);

    template <Element::CODE code> void runElement(SubKernel::Action action, size_t interval);
    template <Element::CODE code> void runBoundary(SubKernel::Action action, size_t region, size_t interval);

    std::vector<HeatTransferElementOperators> subkernels;
    std::vector<std::vector<HeatTransferBoundaryOperators> > boundary;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_H_ */
