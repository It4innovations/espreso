
#ifndef SRC_FETI_FETI_H_
#define SRC_FETI_FETI_H_

#include "config/ecf/linearsolver/feti.h"
#include "esinfo/stepinfo.h"
#include "analysis/pattern/decomposition.feti.h"
#include "math/primitives/vector_dense.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_dnblas.h"
#include "gpu/gpu_spblas.h"

#include <complex>

namespace espreso {

template <typename T> class IterativeSolver;
template <typename T> struct Preconditioner;
template <typename T> struct Projector;
template <typename T> class DualOperator;

template<typename T>
struct FETI {
    struct SystemInfo {
        int domains, clusters;
    };

    FETI(FETIConfiguration &configuration);
    ~FETI();

    void check();

    bool set(const step::Step &step);
    bool update(const step::Step &step);
    bool solve(const step::Step &step);

    FETIConfiguration &configuration;
    SystemInfo sinfo;

    DecompositionFETI *decomposition;
    std::vector<Matrix_CSR<T> > K, assembledK;
    std::vector<Vector_Dense<T> > f, x, BtL, ineqBtL;
    std::vector<Matrix_Dense<T> > MoorePenroseInv;

    std::vector<Matrix_Dense<T> > R1, R2, KR1, KR2;
    std::vector<Matrix_CSR<T> > RegMat;

    std::vector<Matrix_CSR<T> > B1;
    std::vector<std::vector<int> > D2C;
    Vector_Dense<T> c, lb, ub;
    struct Lambdas {
        struct interval { int halo, size; };

        std::vector<interval> intervals;

        // general
        int equalities, size;
        std::vector<int> cmap; // size, ndomains <d0, d1, ..., dn>; size, ndomains <>; ...;
    } lambdas;

    struct Cluster {
        int gl_size = 0;
    } cluster;

    struct {
        bool K = true, B = true;
    } updated;

    IterativeSolver<T> *iterativeSolver = nullptr;
    Preconditioner<T> *preconditioner = nullptr;
    Projector<T> *projector = nullptr;
    DualOperator<T> *dualOperator = nullptr;

    gpu::mgm::device device;
    gpu::mgm::queue main_q;
    std::vector<gpu::mgm::queue> queues;
    std::vector<gpu::dnblas::handle> handles_dense;
    std::vector<gpu::spblas::handle> handles_sparse;
};

}

#endif /* SRC_FETI_FETI_H_ */
