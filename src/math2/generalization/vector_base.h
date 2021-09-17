
#ifndef SRC_MATH2_GENERALIZATION_VECTOR_BASE_H_
#define SRC_MATH2_GENERALIZATION_VECTOR_BASE_H_

#include "analysis/composer/elementmapping.h"
#include "math2/primitives/vector_sparse.h"

#include <complex>
#include <vector>

namespace espreso {

template <typename T> class Vector_Dense;
template <template<typename> typename Vector, typename T> class Vector_Distributed;

template <typename V, typename T> class Vector_Base_Common {
public:
	Vector_Base_Common(): touched(false) {}
	virtual ~Vector_Base_Common() {};

	virtual V* copyPattern() =0;
	virtual void store(const char *file) =0;
	virtual void store(std::vector<double> &output) =0;

	virtual void set(const T &value) =0;
	virtual void scale(const T&alpha) =0;

	virtual void copy(const V *in) =0;
	virtual void add(const T &alpha, const V *a) =0;

	virtual T norm() =0;
	virtual T max() =0;
	virtual T absmax() =0;
	virtual T dot(const V *other) =0;

	virtual void addTo(const T &alpha, Vector_Distributed<Vector_Dense , T> *a) const =0;
	virtual void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, T> *a) const =0;

	ElementMapping<T> mapping;
	bool touched;
};

template <typename T>
class Vector_Base: public Vector_Base_Common<Vector_Base<T>, T>
{
public:
	using Vector_Base_Common<Vector_Base<T>, T>::copy;
	using Vector_Base_Common<Vector_Base<T>, T>::add;
	using Vector_Base_Common<Vector_Base<T>, T>::addTo;

	virtual void copy(const Vector_Base<T> *in, int offset, int size, int step) =0;
	virtual void add(const T &alpha, const Vector_Base<T> *a, int offset, int size, int step) =0;

	virtual void addTo(const T &alpha, Vector_Distributed<Vector_Dense, T> *a, int offset, int size, int step) const =0;
	virtual void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, T> *a, int offset, int size, int step) const =0;
};

template <typename T>
class Vector_Base<std::complex<T> >: public Vector_Base_Common<Vector_Base<std::complex<T> >, std::complex<T> >
{
public:
	virtual void copyReal(const Vector_Base<T> *in) =0;
	virtual void copyImag(const Vector_Base<T> *in) =0;
	virtual void copyRealTo(Vector_Base<T> *in) const =0;
	virtual void copyImagTo(Vector_Base<T> *in) const =0;
};

}

#endif /* SRC_MATH2_GENERALIZATION_VECTOR_BASE_H_ */
