
#ifndef SRC_MATH2_FETI_VECTOR_KERNEL_H_
#define SRC_MATH2_FETI_VECTOR_KERNEL_H_

#include "esinfo/envinfo.h"
#include "math/primitives/vector_dense.h"
#include "wrappers/mpi/communication.h"

#include <vector>
#include <memory>

namespace espreso {

template <typename T>
struct Vector_Kernel: public Vector_Dense<T, int> {

    enum: size_t {
        align = 64U
    };

    static void set(int offset, int size, int totalSize);

    Vector_Kernel();

    void resize();
    void resize(int size);
    void synchronize();

    T dot(const Vector_Dense<T> &other) const;
    T dot() const;

    using Vector_Dense<T>::size;
    using Vector_Dense<T>::vals;

    static int offset, localSize, totalSize;
    static std::vector<size_t> distribution;
};

template <typename T> int                 Vector_Kernel<T>::offset = 0;
template <typename T> int                 Vector_Kernel<T>::localSize = 0;
template <typename T> int                 Vector_Kernel<T>::totalSize = 0;
template <typename T> std::vector<size_t> Vector_Kernel<T>::distribution = { 0, 0 };

}

#endif /* SRC_MATH2_FETI_VECTOR_KERNEL_H_ */
