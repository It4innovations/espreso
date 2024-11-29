
#ifndef SRC_ANALYSIS_MATH_MATRIX_DISTRIBUTED_H_
#define SRC_ANALYSIS_MATH_MATRIX_DISTRIBUTED_H_

#include "analysis/math/math.physics.h"
#include "analysis/math/matrix_base.h"
#include "analysis/math/vector_distributed.h"
#include "analysis/pattern/apply.h"
#include "analysis/pattern/synchronization.h"
#include "analysis/pattern/decomposition.direct.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"

#include <vector>

namespace espreso {

template <typename T>
class Matrix_Distributed: public Matrix_Base<T> {
public:
    void synchronize()
    {
        sync->gatherFromUpper(*static_cast<Matrix_Distributed<T>*>(this));
    }

    Matrix_Base<T>* copyPattern()
    {
        Matrix_Distributed<T> *m = new Matrix_Distributed<T>();
        m->type = m->cluster.type = this->type;
        m->shape = m->cluster.shape = this->shape;
        m->cluster.pattern(cluster);
        m->decomposition = this->decomposition;
        m->applyMatrix = this->applyMatrix;
        m->sync = this->sync;
        return m;
    }

    void store(const char *file)
    {
        math::store(*static_cast<Matrix_Distributed<T>*>(this), file);
    }

    Matrix_Base<T>* set(const T &value)
    {
        math::set(cluster, value);
        return this;
    }

    Matrix_Base<T>* copy(const Matrix_Base<T> *in, const Selection &rows = Selection(), const Selection &cols = Selection())
    {
        in->copyTo(static_cast<Matrix_Distributed<T>*>(this), rows, cols);
        return this;
    }

    Matrix_Base<T>* add(const T &alpha, const Matrix_Base<T> *a, const Selection &rows = Selection(), const Selection &cols = Selection())
    {
        a->addTo(alpha, static_cast<Matrix_Distributed<T>*>(this), rows, cols);
        return this;
    }

    esint size()
    {
        return cluster.nnz;
    }

    void apply(const T &alpha, const Vector_Base<T> *in, const T &beta, Vector_Base<T> *out)
    {
        if (dynamic_cast<const Vector_Distributed<Vector_Dense, T>*>(in) && dynamic_cast<Vector_Distributed<Vector_Dense, T>*>(out)) {
            applyMatrix->apply(*this, *static_cast<Vector_Distributed<Vector_Dense, T>*>(out), alpha, beta, *static_cast<const Vector_Distributed<Vector_Dense, T>*>(in));
        } else {
            eslog::error("call empty function Matrix_Distributed::apply\n");
        }
    }

    void copyTo(Matrix_Distributed<T> *a, const Selection &rows = Selection(), const Selection &cols = Selection()) const
    {
        math::copy(a->cluster, this->cluster, rows, cols);
    }

    void copyTo(Matrix_FETI<T> *a, const Selection &rows = Selection(), const Selection &cols = Selection()) const
    {
        eslog::error("call empty function\n");
    }

    void addTo(const T &alpha, Matrix_Distributed<T> *a, const Selection &rows = Selection(), const Selection &cols = Selection()) const
    {
        math::add(a->cluster, alpha, this->cluster, rows, cols);
    }

    void addTo(const T &alpha, Matrix_FETI<T> *a, const Selection &rows = Selection(), const Selection &cols = Selection()) const
    {
        eslog::error("call empty function\n");
    }

    void printEigenValues(const char *name, size_t n) const
    {
        if (info::mpi::size > 1) {
            eslog::error("cannot print eigen values of matrix distributed to more MPI processes.\n");
        }
        Matrix_Dense<T> m; m.resize(cluster.nrows, cluster.ncols); math::set(m, T{0});
        Vector_Dense<T> e; e.resize(cluster.nrows);
        for (esint r = 0; r < cluster.nrows; ++r) {
            for (esint c = cluster.rows[r]; c < cluster.rows[r + 1]; ++c) {
                m.vals[r * cluster.ncols + cluster.cols[c - Indexing::CSR]] = cluster.vals[c - Indexing::CSR];
            }
        }
        math::lapack::get_eig_sym(m, e);
        printf("%s:", name);
        for (size_t v = 0; v < n; ++v) {
            printf(" %+.5e", e.vals[v]);
        }
        if (n < (size_t)e.size) {
            printf(" ... %+.5e", e.vals[e.size - 1]);
        }
        printf("\n");
    }

    Matrix_CSR<T, esint> cluster;
    DecompositionDirect *decomposition;

    ApplyMatrix<Matrix_Distributed<T>, T > *applyMatrix;
    Synchronization<Matrix_Distributed<T> > *sync;
};

}

#endif /* SRC_ANALYSIS_MATH_MATRIX_DISTRIBUTED_H_ */
