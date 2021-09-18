
#ifndef SRC_MATH2_GENERALIZATION_VECTOR_BASE_H_
#define SRC_MATH2_GENERALIZATION_VECTOR_BASE_H_

#include "analysis/composer/elementmapping.h"

#include <complex>
#include <vector>

namespace espreso {

template <typename T> class Vector_Base;
template <typename T> class Vector_Dense;
template <typename T> class Vector_Sparse;
template <template<typename> typename Vector, typename T> class Vector_Distributed;

template <typename T> class Vector_Base_Common {
public:
	Vector_Base_Common(): touched(false) {}
	virtual ~Vector_Base_Common() {};

	virtual Vector_Base<T>* copyPattern() =0;
	virtual void store(const char *file) =0;
	virtual void store(std::vector<double> &output) =0;

	virtual void set(const T &value) =0;
	virtual void scale(const T&alpha) =0;

	virtual void copy(const Vector_Base<T> *in) =0;
	virtual void add(const T &alpha, const Vector_Base<T> *a) =0;

	virtual T norm() =0;
	virtual T max() =0;
	virtual T absmax() =0;
	virtual T dot(const Vector_Base<T> *other) =0;

	virtual void copyTo(Vector_Distributed<Vector_Dense , T> *a) const =0;
	virtual void copyTo(Vector_Distributed<Vector_Sparse, T> *a) const =0;
	virtual void addTo(const T &alpha, Vector_Distributed<Vector_Dense , T> *a) const =0;
	virtual void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, T> *a) const =0;

	ElementMapping<T> mapping;
	bool touched;
};

template <typename T>
class Vector_Base: public Vector_Base_Common<T>
{
public:
	using Vector_Base_Common<T>::copy;
	using Vector_Base_Common<T>::add;
	using Vector_Base_Common<T>::copyTo;
	using Vector_Base_Common<T>::addTo;

	virtual void copyReal(const Vector_Distributed<Vector_Dense , std::complex<T> > *a) =0;
	virtual void copyReal(const Vector_Distributed<Vector_Sparse, std::complex<T> > *a) =0;
	virtual void copyImag(const Vector_Distributed<Vector_Dense , std::complex<T> > *a) =0;
	virtual void copyImag(const Vector_Distributed<Vector_Sparse, std::complex<T> > *a) =0;

	virtual void copyToReal(Vector_Distributed<Vector_Dense , std::complex<T> > *a) const =0;
	virtual void copyToReal(Vector_Distributed<Vector_Sparse, std::complex<T> > *a) const =0;
	virtual void copyToImag(Vector_Distributed<Vector_Dense , std::complex<T> > *a) const =0;
	virtual void copyToImag(Vector_Distributed<Vector_Sparse, std::complex<T> > *a) const =0;

	virtual void copy(const Vector_Base<T> *in, int offset, int size, int step) =0;
	virtual void add(const T &alpha, const Vector_Base<T> *a, int offset, int size, int step) =0;

	virtual void copyTo(Vector_Distributed<Vector_Dense, T> *a, int offset, int size, int step) const =0;
	virtual void copyTo(Vector_Distributed<Vector_Sparse, T> *a, int offset, int size, int step) const =0;
	virtual void addTo(const T &alpha, Vector_Distributed<Vector_Dense, T> *a, int offset, int size, int step) const =0;
	virtual void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, T> *a, int offset, int size, int step) const =0;
};

template <typename T>
class Vector_Base<std::complex<T> >: public Vector_Base_Common<std::complex<T> >
{
public:
	virtual void copyReal(const Vector_Base<T> *in) =0;
	virtual void copyImag(const Vector_Base<T> *in) =0;
	virtual void copyRealTo(Vector_Base<T> *in) const =0;
	virtual void copyImagTo(Vector_Base<T> *in) const =0;
};

}

#endif /* SRC_MATH2_GENERALIZATION_VECTOR_BASE_H_ */
