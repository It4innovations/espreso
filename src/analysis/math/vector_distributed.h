
#ifndef SRC_ANALYSIS_MATH_VECTOR_DISTRIBUTED_H_
#define SRC_ANALYSIS_MATH_VECTOR_DISTRIBUTED_H_

#include "analysis/math/math.physics.h"
#include "analysis/math/vector_base.h"
#include "analysis/math/vector_feti.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "analysis/pattern/decomposition.direct.h"
#include "analysis/pattern/synchronization.h"
#include "wrappers/mpi/communication.h"

#include <vector>

namespace espreso {

template <template<typename, typename, typename> typename Vector, typename T>
class Vector_Distributed: public Vector_Base<T> {
public:
    void synchronize()
    {
        sync->synchronize(*static_cast<Vector_Distributed<Vector, T>*>(this));
    }

    void scatter()
    {
        sync->scatterToUpper(*static_cast<Vector_Distributed<Vector, T>*>(this));
    }

    Vector_Distributed<Vector, T>* copyPattern()
    {
        Vector_Distributed<Vector, T> *v = new Vector_Distributed<Vector, T>();
        v->cluster.pattern(cluster);
        v->decomposition = this->decomposition;
        v->sync = this->sync;
        return v;
    }

    void store(const char *file)
    {
        math::store(*static_cast<Vector_Distributed<Vector, T>*>(this), file);
    }

    void storeTo(std::vector<double> &output)
    {
        for (size_t i = 0; i < output.size(); ++i) {
            output[i] = this->cluster.vals[i];
        }
    }

    void setFrom(std::vector<double> &output)
    {
        for (size_t i = 0; i < output.size(); ++i) {
            this->cluster.vals[i] = output[i];
        }
    }

    Vector_Base<T>* set(const T &value)
    {
        math::set(cluster, value);
        return this;
    }

    Vector_Base<T>* scale(const T &value)
    {
        math::scale(value, cluster);
        return this;
    }

    Vector_Base<T>* copy(const Vector_Base<T> *a, const Selection &rows = Selection())
    {
        a->copyTo(static_cast<Vector_Distributed<Vector, T>*>(this), rows);
        return this;
    }

    Vector_Base<T>* add(const T &alpha, const Vector_Base<T> *a, const Selection &rows = Selection())
    {
        a->addTo(alpha, static_cast<Vector_Distributed<Vector, T>*>(this));
        return this;
    }

    T norm()
    {
        T dot = math::blas::dot(cluster.size - decomposition->halo.size(), cluster.vals + decomposition->halo.size(), 1, cluster.vals + decomposition->halo.size(), 1);
        Communication::allReduce(&dot, NULL, 1, MPI_DOUBLE, MPI_SUM);
        return std::sqrt(dot);
    }

    T max()
    {
        eslog::error("call empty function: max\n");
        return 0;
    }

    T absmax()
    {
        eslog::error("call empty function: absmax\n");
        return 0;
    }

    T dot(const Vector_Base<T> *other)
    {
        eslog::error("call empty function: dot\n");
        return 0;
    }

    void copyTo(Vector_Distributed<Vector_Dense , T> *a, const Selection &rows = Selection()) const
    {
        math::copy(a->cluster, cluster, rows);
    }

    void copyTo(Vector_Distributed<Vector_Sparse, T> *a, const Selection &rows = Selection()) const
    {
        math::copy(a->cluster, cluster, rows);
    }

    void copyTo(Vector_FETI<Vector_Dense , T> *a, const Selection &rows = Selection()) const
    {
        esint index = 0;
        for (auto dmap = a->decomposition->dmap->cbegin(); dmap != a->decomposition->dmap->cend(); ++dmap, ++index) {
            for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                if (a->decomposition->ismy(di->domain)) {
                    a->domains[di->domain - a->decomposition->dbegin].vals[di->index] = cluster.vals[index];
                }
            }
        }
    }

    void spliteTo(Vector_FETI<Vector_Dense , T> *a, const Selection &rows = Selection()) const
    {
        esint index = 0;
        for (auto dmap = a->decomposition->dmap->cbegin(); dmap != a->decomposition->dmap->cend(); ++dmap, ++index) {
            for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                if (a->decomposition->ismy(di->domain)) {
                    a->domains[di->domain - a->decomposition->dbegin].vals[di->index] = cluster.vals[index] / dmap->size();
                }
            }
        }
    }

    void copyTo(Vector_FETI<Vector_Sparse, T> *a, const Selection &rows = Selection()) const
    {
        eslog::error("call empty function\n");
    }

    void addTo(const T &alpha, Vector_Distributed<Vector_Dense, T> *a, const Selection &rows = Selection()) const
    {
        math::add(a->cluster, alpha, cluster, rows);
    }

    void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, T> *a, const Selection &rows = Selection()) const
    {
        math::add(a->cluster, alpha, cluster, rows);
    }

    void addTo(const T &alpha, Vector_FETI<Vector_Dense, T> *a, const Selection &rows = Selection()) const
    {
        eslog::error("call empty function\n");
    }

    void addTo(const T &alpha, Vector_FETI<Vector_Sparse, T> *a, const Selection &rows = Selection()) const
    {
        eslog::error("call empty function\n");
    }

    void print() const
    {
        for (esint n = 0; n < info::mesh->nodes->size; ++n) {
            printf("%2d [%+.14e %+.14e %+.14e]\n", n, cluster.vals[3 * n + 0], cluster.vals[3 * n + 1], cluster.vals[3 * n + 2]);
        }
    }

    Vector<T, esint, cpu_allocator> cluster;
    DecompositionDirect *decomposition;
    Synchronization<Vector_Distributed<Vector, T> > *sync;
};

}

#endif /* SRC_ANALYSIS_MATH_VECTOR_DISTRIBUTED_H_ */
