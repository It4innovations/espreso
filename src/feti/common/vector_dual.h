
#ifndef SRC_MATH2_FETI_VECTOR_DUAL_H_
#define SRC_MATH2_FETI_VECTOR_DUAL_H_

#include "math/primitives/vector_dense.h"
#include "feti/feti.h"
#include "dual_buffer.h"

#include <vector>

namespace espreso {

struct FETIDecomposition;

template <typename T>
struct Vector_Dual: public Vector_Dense<T> {

	template <typename Type> friend struct Matrix_Dual_Orthogonal;

	void resize(int nhalo, int size);
	void synchronize(Dual_Buffer<T> &buffer);
	void copyToWithoutHalo(Vector_Dense<T> &to);
	T dot(const Vector_Dense<T> &other) const;
	T dot() const;

	int nhalo;
	using Vector_Dense<T>::size;
	using Vector_Dense<T>::vals;
};

}

#endif /* SRC_MATH2_FETI_VECTOR_DUAL_H_ */
