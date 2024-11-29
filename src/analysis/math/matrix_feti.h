
#ifndef SRC_ANALYSIS_MATH_MATRIX_FETI_H_
#define SRC_ANALYSIS_MATH_MATRIX_FETI_H_

#include "analysis/math/math.physics.h"
#include "analysis/math/matrix_base.h"
#include "analysis/math/vector_feti.h"
#include "analysis/math/vector_distributed.h"
#include "analysis/pattern/apply.h"
#include "analysis/pattern/decomposition.feti.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/math.h"

#include <vector>

namespace espreso {

template <typename T>
class Matrix_FETI: public Matrix_Base<T> {
public:
    void synchronize()
    {

    }

    Matrix_Base<T>* copyPattern()
    {
        Matrix_FETI<T> *m = new Matrix_FETI<T>();
        m->type = this->type;
        m->shape = this->shape;
        m->domains.resize(domains.size());
        #pragma omp parallel for
        for (size_t d = 0; d < domains.size(); ++d) {
            m->domains[d].type = domains[d].type;
            m->domains[d].shape = domains[d].shape;
            m->domains[d].pattern(domains[d]);
        }
        m->decomposition = decomposition;
        m->applyMatrix = applyMatrix;
        return m;
    }

    void store(const char *file)
    {
        math::store(*static_cast<Matrix_FETI<T>*>(this), file);
    }

    Matrix_Base<T>* set(const T &value)
    {
        #pragma omp parallel for
        for (size_t d = 0; d < this->domains.size(); ++d) {
            math::set(this->domains[d], value);
        }
        return this;
    }

    Matrix_Base<T>* copy(const Matrix_Base<T> *in, const Selection &rows = Selection(), const Selection &cols = Selection())
    {
        in->copyTo(static_cast<Matrix_FETI<T>*>(this), rows, cols);
        return this;
    }

    Matrix_Base<T>* add(const T &alpha, const Matrix_Base<T> *a, const Selection &rows = Selection(), const Selection &cols = Selection())
    {
        a->addTo(alpha, static_cast<Matrix_FETI<T>*>(this), rows, cols);
        return this;
    }

    esint size()
    {
        return 0; // how to compute size?
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
        eslog::error("call empty function\n");
    }

    void copyTo(Matrix_FETI<T> *a, const Selection &rows = Selection(), const Selection &cols = Selection()) const
    {
        #pragma omp parallel for
        for (size_t d = 0; d < this->domains.size(); ++d) {
            math::copy(a->domains[d], this->domains[d], rows, cols);
        }
    }

    void addTo(const T &alpha, Matrix_Distributed<T> *a, const Selection &rows = Selection(), const Selection &cols = Selection()) const
    {
        eslog::error("call empty function\n");
    }

    void addTo(const T &alpha, Matrix_FETI<T> *a, const Selection &rows = Selection(), const Selection &cols = Selection()) const
    {
        #pragma omp parallel for
        for (size_t d = 0; d < this->domains.size(); ++d) {
            math::add(a->domains[d], alpha, this->domains[d], rows, cols);
        }
    }

    void printEigenValues(const char *name, size_t n) const
    {
        std::vector<Vector_Dense<T> > e(domains.size());
        for (size_t d = 0; d < domains.size(); ++d) {
            Matrix_Dense<T> m; m.resize(domains[d].nrows, domains[d].ncols); math::set(m, T{0});
            e[d].resize(domains[d].nrows);

            for (esint r = 0; r < domains[d].nrows; ++r) {
                for (esint c = domains[d].rows[r]; c < domains[d].rows[r + 1]; ++c) {
                    m.vals[r * domains[d].ncols + domains[d].cols[c - Indexing::CSR]] = domains[d].vals[c - Indexing::CSR];
                }
            }
            math::lapack::get_eig_sym(m, e[d]);
        }
        Communication::serialize([&] () {
            for (size_t d = 0; d < domains.size(); ++d) {
                printf("%s:", name);
                for (size_t v = 0; v < n; ++v) {
                    printf(" %+.5e", e[d].vals[v]);
                }
                if (n < (size_t)e[d].size) {
                    printf(" ... %+.5e", e[d].vals[e[d].size - 1]);
                }
                printf("\n");
            }
        });
    }

    std::vector<Matrix_CSR<T, int> > domains;
    DecompositionFETI *decomposition;

    ApplyMatrix<Matrix_FETI<T>, T> *applyMatrix;
};

}

#endif /* SRC_ANALYSIS_MATH_MATRIX_FETI_H_ */
