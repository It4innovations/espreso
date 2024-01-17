
#ifndef SRC_MATH2_FETI_VECTOR_DUAL_H_
#define SRC_MATH2_FETI_VECTOR_DUAL_H_

#include "math/primitives/vector_dense.h"

#include <vector>

namespace espreso {

struct DOFsDecomposition;

template <typename T>
struct Vector_Dual: public Vector_Dense<T, int> {

	template <typename Type> friend struct Matrix_Dual_Orthogonal;

	static void set(esint dirichlet, esint nhalo, const std::vector<esint> &cmap, const DOFsDecomposition &decomposition);

	void resize();
	void synchronize();
	void copyTo(Vector_Dense<T, int> &to) const;
	void copyToWithoutHalo(Vector_Dense<T, int> &to) const;
	void scale(const T &alpha);
	void add(const T &alpha, const Vector_Dense<T, int> &other);

	double dot(const Vector_Dense<T, int> &other) const;
	double dot() const;

protected:
	static esint dirichlet, nhalo, size;
	static std::vector<int> nmap, neighbors;
	static std::vector<std::vector<T> > sBuffer, rBuffer;

};

template <typename T> esint Vector_Dual<T>::dirichlet;
template <typename T> esint Vector_Dual<T>::nhalo;
template <typename T> esint Vector_Dual<T>::size;
template <typename T> std::vector<int> Vector_Dual<T>::nmap;
template <typename T> std::vector<int> Vector_Dual<T>::neighbors;
template <typename T> std::vector<std::vector<T> > Vector_Dual<T>::sBuffer;
template <typename T> std::vector<std::vector<T> > Vector_Dual<T>::rBuffer;

}

#endif /* SRC_MATH2_FETI_VECTOR_DUAL_H_ */
