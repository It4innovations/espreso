
#ifndef SRC_MATH2_FETI_VECTOR_DUAL_H_
#define SRC_MATH2_FETI_VECTOR_DUAL_H_

#include "math/primitives/vector_dense.h"
#include "math/feti/lmap.h"

#include <vector>

namespace espreso {

template <typename T>
struct Vector_Dual: public Vector_Dense<T> {

	template <typename Type> friend struct Matrix_Dual_Orthogonal;

	static void set(esint nhalo, esint size, const std::vector<LMAP> &lmap, const std::vector<int> &neighbors);

	void resize();
	void synchronize();
	void copyTo(Vector_Dense<T> &to) const;
	void copyToWithoutHalo(Vector_Dense<T> &to) const;
	void scale(const T &alpha);
	void add(const T &alpha, const Vector_Dense<T> &other);

	double dot(const Vector_Dense<T> &other) const;
	double dot() const;

protected:
	static esint nhalo, nshared, localSize;
	static Vector_Dense<T> halo;
	static std::vector<esint> distribution;
	static std::vector<int> neighbors;
};

template <typename T>
esint Vector_Dual<T>::nhalo = 0;
template <typename T>
esint Vector_Dual<T>::nshared = 0;
template <typename T>
esint Vector_Dual<T>::localSize = 0;
template <typename T>
Vector_Dense<T> Vector_Dual<T>::halo;
template <typename T>
std::vector<esint> Vector_Dual<T>::distribution;
template <typename T>
std::vector<esint> Vector_Dual<T>::neighbors;

}

#endif /* SRC_MATH2_FETI_VECTOR_DUAL_H_ */
