
#ifndef SRC_ANALYSIS_MATH_VECTOR_FETI_H_
#define SRC_ANALYSIS_MATH_VECTOR_FETI_H_

#include "analysis/math/math.physics.h"
#include "analysis/math/vector_base.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "analysis/pattern/decomposition.feti.h"
#include "wrappers/mpi/communication.h"

#include <vector>

namespace espreso {

template <template<typename, typename, typename> typename Vector, typename T>
class Vector_FETI: public Vector_Base<T> {

    void _store(double &out, const double &in) { out = in; }
    void _store(double &out, const std::complex<double> &in) { out = in.real(); }

public:
    void synchronize()
    {
//        _sync->synchronize(*static_cast<Vector_FETI<Vector, T>*>(this));
    }

    void split()
    {
//        _sync->split(*static_cast<Vector_FETI<Vector, T>*>(this));
    }

    void scatter()
    {
//        synchronization->scatterToUpper(*static_cast<Vector_FETI<Vector, T>*>(this));
    }

    Vector_FETI<Vector, T>* copyPattern()
    {
        Vector_FETI<Vector, T> *m = new Vector_FETI<Vector, T>();
        m->domains.resize(domains.size());
        #pragma omp parallel for
        for (size_t d = 0; d < domains.size(); ++d) {
            m->domains[d].pattern(domains[d]);
        }
        m->decomposition = decomposition;
//        m->_sync = _sync;
        return m;
    }

    void store(const char *file)
    {
        math::store(*static_cast<Vector_FETI<Vector, T>*>(this), file);
    }

    void storeTo(std::vector<double> &output)
    {
        auto dmap = decomposition->dmap->cbegin();
        for (size_t i = 0; i < output.size(); ++i, ++dmap) {
            for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                if (decomposition->ismy(di->domain)) {
                    _store(output[i], domains[di->domain - decomposition->dbegin].vals[di->index]);
                    break; // we assume synchronization inside the solver
                }
            }
        }
    }

    void setFrom(std::vector<double> &output)
    {
        auto dmap = decomposition->dmap->cbegin();
        for (size_t i = 0; i < output.size(); ++i, ++dmap) {
            for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                if (decomposition->ismy(di->domain)) {
                    _store(domains[di->domain - decomposition->dbegin].vals[di->index], output[i]);
                }
            }
        }
    }

    Vector_Base<T>* set(const T &value)
    {
        #pragma omp parallel for
        for (size_t d = 0; d < domains.size(); ++d) {
            math::set(domains[d], value);
        }
        return this;
    }

    Vector_Base<T>* scale(const T &value)
    {
        #pragma omp parallel for
        for (size_t d = 0; d < domains.size(); ++d) {
            math::scale(value, domains[d]);
        }
        return this;
    }

    Vector_Base<T>* copy(const Vector_Base<T> *in, const Selection &rows = Selection())
    {
        in->copyTo(static_cast<Vector_FETI<Vector, T>*>(this), rows);
        return this;
    }

    Vector_Base<T>* add(const T &alpha, const Vector_Base<T> *a, const Selection &rows = Selection())
    {
        a->addTo(alpha, static_cast<Vector_FETI<Vector, T>*>(this), rows);
        return this;
    }

    esint size()
    {
        return decomposition->halo.size() + (decomposition->end - decomposition->begin);
    }

    T _dot(const std::vector<Vector_Sparse<T> > &domains)
    {
        eslog::error("call empty function: max\n");
        return T{};
    }

    T _dot(const std::vector<Vector_Dense<T> > &domains)
    {
        T dot = 0;
        for (auto map = decomposition->dmap->cbegin(); map != decomposition->dmap->cend(); ++map) {
            if (decomposition->ismy(map->begin()->domain)) {
                const T &val = domains[map->begin()->domain - decomposition->dbegin].vals[map->begin()->index];
                dot += val * val;
            }
        }
        return dot;
    }

    T norm()
    {
        T dot = _dot(domains);
        Communication::allReduce(&dot, NULL, 1, MPITools::getType(dot).mpitype, MPI_SUM);
        return std::sqrt(dot);
    }

    T max()
    {
        eslog::error("call empty function: max\n");
        return T{};
    }

    T absmax()
    {
        eslog::error("call empty function: absmax\n");
        return T{};
    }

    T dot(const Vector_Base<T> *other)
    {
        eslog::error("call empty function: dot\n");
        return T{};
    }

    void copyTo(Vector_Distributed<Vector_Dense , T> *a, const Selection &rows = Selection()) const
    {
        esint index = 0;
        for (auto dmap = decomposition->dmap->cbegin(); dmap != decomposition->dmap->cend(); ++dmap, ++index) {
            for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                if (decomposition->ismy(di->domain)) {
                    a->cluster.vals[index] = domains[di->domain - decomposition->dbegin].vals[di->index];
                    break;
                }
            }
        }
    }

    void copyTo(Vector_Distributed<Matrix_Dense , T> *a, const Selection &rows = Selection()) const
    {
        esint index = 0;
        for (auto dmap = decomposition->dmap->cbegin(); dmap != decomposition->dmap->cend(); ++dmap, ++index) {
            for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                if (decomposition->ismy(di->domain)) {
                    for (int c = 0; c < a->cluster.ncols; ++c) {
                        a->cluster.vals[index + c * a->cluster.nrows] = domains[di->domain - decomposition->dbegin].vals[di->index + c * a->cluster.nrows];
                    }
                    break;
                }
            }
        }
    }

    void averageTo(Vector_Distributed<Vector_Dense , T> *a, const Selection &rows = Selection()) const
    {
        sumTo(a, rows);
        esint index = 0;
        for (auto dmap = decomposition->dmap->cbegin(); dmap != decomposition->dmap->cend(); ++dmap, ++index) {
            a->cluster.vals[index] /= dmap->size();
        }
    }

    void sumTo(Vector_Distributed<Vector_Dense , T> *a, const Selection &rows = Selection()) const
    {
        math::set(a->cluster, 0.);
        esint index = 0;
        for (auto dmap = decomposition->dmap->cbegin(); dmap != decomposition->dmap->cend(); ++dmap, ++index) {
            for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                if (decomposition->ismy(di->domain)) {
                    a->cluster.vals[index] += domains[di->domain - decomposition->dbegin].vals[di->index];
                }
            }
        }
        a->synchronize();
    }

    void copyTo(Vector_Distributed<Vector_Sparse, T> *a, const Selection &rows = Selection()) const
    {
        eslog::error("call empty function\n");
    }

    void copyTo(Vector_FETI<Vector_Dense , T> *a, const Selection &rows = Selection()) const
    {
        #pragma omp parallel for
        for (size_t d = 0; d < domains.size(); ++d) {
            math::copy(a->domains[d], domains[d], rows);
        }
    }

    void copyTo(Vector_FETI<Matrix_Dense , T> *a, const Selection &rows = Selection()) const
    {
        #pragma omp parallel for
        for (size_t d = 0; d < domains.size(); ++d) {
            math::copy(a->domains[d], domains[d], rows);
        }
    }

    void copyTo(Vector_FETI<Vector_Sparse, T> *a, const Selection &rows = Selection()) const
    {
        #pragma omp parallel for
        for (size_t d = 0; d < domains.size(); ++d) {
            math::copy(a->domains[d], domains[d], rows);
        }
    }

    void addTo(const T &alpha, Vector_Distributed<Vector_Dense, T> *a, const Selection &rows = Selection()) const
    {
        auto dmap = decomposition->dmap->cbegin();
        for (esint i = 0; i < a->cluster.size; ++i, ++dmap) {
            for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                if (decomposition->ismy(di->domain)) {
                    a->cluster.vals[i] += alpha * domains[di->domain - decomposition->dbegin].vals[di->index];
                    break;
                }
            }
        }
    }

    void addTo(const T &alpha, Vector_Distributed<Matrix_Dense, T> *a, const Selection &rows = Selection()) const
    {
        auto dmap = decomposition->dmap->cbegin();
        for (esint i = 0; i < a->cluster.nrows; ++i, ++dmap) {
            for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                if (decomposition->ismy(di->domain)) {
                    for (int c = 0; c < a->cluster.ncols; ++c) {
                        a->cluster.vals[i + c * a->cluster.nrows] += alpha * domains[di->domain - decomposition->dbegin].vals[di->index + c * a->cluster.nrows];
                    }
                    break;
                }
            }
        }
    }

    void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, T> *a, const Selection &rows = Selection()) const
    {
        auto dmap = decomposition->dmap->cbegin();
        for (esint i = 0, prev = 0; i < a->cluster.nnz; prev = a->cluster.indices[i++]) {
            dmap += a->cluster.indices[i] - prev;
            for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                if (decomposition->ismy(di->domain)) {
                    a->cluster.vals[i] += alpha * domains[di->domain - decomposition->dbegin].vals[di->index];
                    break;
                }
            }
        }
    }

    void addTo(const T &alpha, Vector_FETI<Vector_Dense, T> *a, const Selection &rows = Selection()) const
    {
        #pragma omp parallel for
        for (size_t d = 0; d < domains.size(); ++d) {
            math::add(a->domains[d], alpha, domains[d], rows);
        }
    }

    void addTo(const T &alpha, Vector_FETI<Matrix_Dense, T> *a, const Selection &rows = Selection()) const
    {
        #pragma omp parallel for
        for (size_t d = 0; d < domains.size(); ++d) {
            math::add(a->domains[d], alpha, domains[d], rows);
        }
    }

    void addTo(const T &alpha, Vector_FETI<Vector_Sparse, T> *a, const Selection &rows = Selection()) const
    {
        #pragma omp parallel for
        for (size_t d = 0; d < domains.size(); ++d) {
            math::add(a->domains[d], alpha, domains[d], rows);
        }
    }

    std::vector<Vector<T, int, cpu_allocator> > domains;
    DecompositionFETI *decomposition;
//    Vector_Sync<Vector, T> *_sync;
};

}

#endif /* SRC_ANALYSIS_MATH_VECTOR_FETI_H_ */
