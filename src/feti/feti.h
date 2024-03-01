
#ifndef SRC_FETI_FETI_H_
#define SRC_FETI_FETI_H_

#include "config/ecf/linearsolver/feti.h"
#include "esinfo/stepinfo.h"
#include "analysis/builder/feti.decomposition.h"
#include "math/primitives/vector_dense.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"

#include <complex>

namespace espreso {

template <typename T> class IterativeSolver;
template <typename T> struct Preconditioner;
template <typename T> struct Projector;
template <typename T> class DualOperator;

template<typename T>
struct FETI {
    struct SystemInfo {
        esint domains, clusters;
        esint R1offset, R2offset;
        esint R1size, R2size;
        esint R1totalSize, R2totalSize;
        size_t lambdasOffset, lambdasLocal, lambdasTotal;
    };

    FETI(FETIConfiguration &configuration);
    ~FETI();

    void info() const;

    bool set(const step::Step &step);
    bool update(const step::Step &step);
    bool solve(const step::Step &step);

    FETIConfiguration &configuration;
    SystemInfo sinfo;

    FETIDecomposition *decomposition;
    std::vector<Matrix_CSR<T> > K;
    std::vector<Vector_Dense<T> > f, x;

    std::vector<Matrix_Dense<T> > R1, R2;
    std::vector<Matrix_CSR<T> > RegMat;

    std::vector<Matrix_CSR<T> > B1;
    std::vector<std::vector<int> > D2C;
    Vector_Dense<T> c, lb, ub;
    struct Lambdas {
        // equality (local)
        esint eq_halo, eq_size;
//        // normal contact
//        esint eq_halo, equalities, size;

        esint equalities, size;
        std::vector<esint> cmap; // size, ndomains <d0, d1, ..., dn>; size, ndomains <>; ...;
    } lambdas;

    struct {
        bool K = true, B = true;
    } updated;

    IterativeSolver<T> *iterativeSolver = nullptr;
    Preconditioner<T> *preconditioner = nullptr;
    Projector<T> *projector = nullptr;
    DualOperator<T> *dualOperator = nullptr;
};

}

#endif /* SRC_FETI_FETI_H_ */
