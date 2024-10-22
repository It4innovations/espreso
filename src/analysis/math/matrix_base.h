
#ifndef SRC_ANALYSIS_MATH_MATRIX_BASE_H_
#define SRC_ANALYSIS_MATH_MATRIX_BASE_H_

#include "selection.h"
#include "vector_base.h"
#include "mapping.h"
#include "math/primitives/matrix_info.h"
#include "math/primitives/matrix_csr.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_ijv.h"

#include <complex>

namespace espreso {

template <typename T> class Matrix_Distributed;
template <typename T> class Matrix_FETI;

template <typename T>
class Matrix_Base {
public:
    Matrix_Base(): type(Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC), shape(Matrix_Shape::FULL), constant(false), filled(false), updated(false) {}

    virtual ~Matrix_Base() {};

    virtual void synchronize() =0;

    virtual Matrix_Base<T>* copyPattern() =0;
    virtual void store(const char *file) =0;

    virtual Matrix_Base<T>* set(const T &value) =0;

    virtual Matrix_Base<T>* copy(const Matrix_Base<T> *in, const Selection &rows = Selection(), const Selection &cols = Selection()) =0;
    virtual void copyTo(Matrix_Distributed<T> *a, const Selection &rows = Selection(), const Selection &cols = Selection()) const =0;
    virtual void copyTo(Matrix_FETI<T> *a, const Selection &rows = Selection(), const Selection &cols = Selection()) const =0;

    virtual Matrix_Base<T>* add(const T &alpha, const Matrix_Base<T> *a, const Selection &rows = Selection(), const Selection &cols = Selection()) =0;
    virtual void addTo(const T &alpha, Matrix_Distributed<T> *a, const Selection &rows = Selection(), const Selection &cols = Selection()) const =0;
    virtual void addTo(const T &alpha, Matrix_FETI<T> *a, const Selection &rows = Selection(), const Selection &cols = Selection()) const =0;

    virtual void apply(const T &alpha, const Vector_Base<T> *in, const T &beta, Vector_Base<T> *out) =0;

    Matrix_Type type;
    Matrix_Shape shape;
    Mapping<T> mapping;
    bool constant, filled, updated;
};

}



#endif /* SRC_ANALYSIS_MATH_MATRIX_BASE_H_ */
