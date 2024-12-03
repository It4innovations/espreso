
#ifndef SRC_FETI_COMMON_MATRIX_DUAL_H_
#define SRC_FETI_COMMON_MATRIX_DUAL_H_

#include "dual_map.h"
#include "feti/feti.h"
#include "math/primitives/matrix_dense.h"

#include <vector>

namespace espreso {

template <typename T>
struct Matrix_Dual: public Matrix_Dense<T> {

    enum: size_t {
        initial_space = 32
    };

    Matrix_Dual()
    {

    }

    Matrix_Dual(int nrhs)
    {
        resize(nrhs);
    }

    virtual ~Matrix_Dual() {}

    virtual void resize()
    {
        // align matrix values ??
        Matrix_Dense<T>::resize(initial_space, Dual_Map::size);
        Matrix_Dense<T>::nrows = 0;
    }

    virtual void resize(int nrhs)
    {
        Matrix_Dense<T>::resize(nrhs, Dual_Map::size);
    }

    virtual void synchronize();

    using Matrix_Dense<T>::nrows;
    using Matrix_Dense<T>::ncols;
    using Matrix_Dense<T>::vals;

    static std::vector<std::vector<T> > send, recv; // buffer
};

template <typename T> std::vector<std::vector<T> > Matrix_Dual<T>::send;
template <typename T> std::vector<std::vector<T> > Matrix_Dual<T>::recv;

}

#endif /* SRC_FETI_COMMON_MATRIX_DUAL_H_ */
