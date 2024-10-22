
#ifndef SRC_ANALYSIS_MATH_VECTOR_BASE_H_
#define SRC_ANALYSIS_MATH_VECTOR_BASE_H_

#include "selection.h"
#include "mapping.h"
#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"

#include <complex>
#include <vector>

namespace espreso {

template <template<typename, typename, typename> typename Vector, typename T> class Vector_Distributed;
template <template<typename, typename, typename> typename Vector, typename T> class Vector_FETI;

template <typename T> class Vector_Base {
public:
    Vector_Base(): constant(false), filled(false), updated(false) {}
    virtual ~Vector_Base() {};

    virtual void synchronize() =0;

    virtual Vector_Base<T>* copyPattern() =0;
    virtual void store(const char *file) =0;
    virtual void storeTo(std::vector<double> &output) =0;
    virtual void setFrom(std::vector<double> &output) =0;

    virtual Vector_Base<T>* set(const T &value) =0;
    virtual Vector_Base<T>* scale(const T &value) =0;

    virtual Vector_Base<T>* copy(const Vector_Base<T> *a, const Selection &rows = Selection()) =0;
    virtual Vector_Base<T>* add(const T &alpha, const Vector_Base<T> *a, const Selection &rows = Selection()) =0;

    virtual T norm() =0;
    virtual T max() =0;
    virtual T absmax() =0;
    virtual T dot(const Vector_Base<T> *other) =0;

    virtual void copyTo(Vector_Distributed<Vector_Dense , T> *a, const Selection &rows = Selection()) const =0;
    virtual void copyTo(Vector_Distributed<Vector_Sparse, T> *a, const Selection &rows = Selection()) const =0;
    virtual void copyTo(Vector_FETI<Vector_Dense , T> *a, const Selection &rows = Selection()) const =0;
    virtual void copyTo(Vector_FETI<Vector_Sparse, T> *a, const Selection &rows = Selection()) const =0;

    virtual void addTo(const T &alpha, Vector_Distributed<Vector_Dense , T> *a, const Selection &rows = Selection()) const =0;
    virtual void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, T> *a, const Selection &rows = Selection()) const =0;
    virtual void addTo(const T &alpha, Vector_FETI<Vector_Dense , T> *a, const Selection &rows = Selection()) const =0;
    virtual void addTo(const T &alpha, Vector_FETI<Vector_Sparse, T> *a, const Selection &rows = Selection()) const =0;

    virtual void print() const =0;

    Mapping<T> mapping;
    bool constant, filled, updated;
};

}

#endif /* SRC_ANALYSIS_MATH_VECTOR_BASE_H_ */
