
#ifndef SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_H_
#define SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_H_

#include "analysis/math/matrix_base.h"
#include "analysis/math/vector_base.h"
#include "analysis/math/vector_feti.h"
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

namespace espreso {

struct HeatTransferConfiguration;
struct HeatTransferLoadStepConfiguration;
struct SteadyState;

class HeatTransfer: public Assembler
{
public:
    HeatTransfer(HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration);

    bool analyze(const step::Step &step);
    void getInitialTemperature(Vector_Base<double> *x);

    void connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet);
    void evaluate(const step::Step &step, step::Time &time, Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet);
    void updateSolution(const step::Step &step, Vector_Distributed<Vector_Dense, double> *x);

    HeatTransferConfiguration &settings;
    HeatTransferLoadStepConfiguration &configuration;

    struct Results {
        static NodeData *temperature, *initialTemperature, *thickness;
        static ElementData *translationMotion, *gradient, *flux;
    };

protected:
    void elements(SubKernel::Action action, const step::Step &step, size_t interval);
    void boundary(SubKernel::Action action, const step::Step &step, size_t region, size_t interval);
    void nodes(SubKernel::Action action, const step::Step &step, size_t region, size_t interval);
    void bem(SubKernel::Action action, size_t domain, double *BETI);

    std::vector<HeatTransferElementOperators> elementKernels;
    std::vector<std::vector<HeatTransferBoundaryOperators> > faceKernels;
    std::vector<std::vector<HeatTransferNodeOperators> > nodeKernels;

    Vector_FETI<Vector_Dense, double> xBEM;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_H_ */
