
#ifndef SRC_FETI_COMMON_VECTOR_DUAL_H_
#define SRC_FETI_COMMON_VECTOR_DUAL_H_

#include "dual_map.h"
#include "feti/feti.h"
#include "math/primitives/vector_dense.h"

#include <vector>

namespace espreso {

template <typename T>
struct Vector_Dual: public Vector_Dense<T> {

    template <typename Type> friend struct Matrix_Dual_Orthogonal;

    Vector_Dual();

    void resize();
    void synchronize();
    void copyToWithoutHalo(Vector_Dense<T> &to) const;
    T dot(const Vector_Dense<T> &other) const;
    T dot() const;

    using Vector_Dense<T>::size;
    using Vector_Dense<T>::vals;

    static void initBuffers()
    {
        send.resize(Dual_Map::nsize.size());
        recv.resize(Dual_Map::nsize.size());
        for (size_t i = 0; i < Dual_Map::nsize.size(); ++i) {
            send[i].resize(Dual_Map::nsize[i]);
            recv[i].resize(Dual_Map::nsize[i]);
        }
    }

    static std::vector<std::vector<T> > send, recv; // buffer
};

template <typename T> std::vector<std::vector<T> > Vector_Dual<T>::send;
template <typename T> std::vector<std::vector<T> > Vector_Dual<T>::recv;

}

#endif /* SRC_FETI_COMMON_VECTOR_DUAL_H_ */
