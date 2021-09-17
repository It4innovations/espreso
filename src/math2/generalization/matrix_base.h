
#ifndef SRC_MATH2_GENERALIZATION_MATRIX_BASE_H_
#define SRC_MATH2_GENERALIZATION_MATRIX_BASE_H_

#include "vector_base.h"

#include "analysis/composer/elementmapping.h"

#include <complex>

namespace espreso {

template <typename M, typename T>
class Matrix_Base_Common {
public:
	Matrix_Base_Common(): type(Matrix_Type::REAL_UNSYMMETRIC), shape(Matrix_Shape::FULL), touched(false) {}

	virtual ~Matrix_Base_Common() {};

	virtual void commit() =0;

//	virtual Matrix_Base* copy() =0;
	virtual M* copyPattern() =0;
	virtual void store(const char *file) =0;

	virtual void set(const T &value) =0;
	virtual void scale(const T &alpha) =0;

	virtual void copy(const M *in) =0;
	virtual void add(const T &alpha, const M *a) =0;
	virtual void apply(const T &alpha, const Vector_Base<T> *in, const T &beta, Vector_Base<T> *out) =0;

	Matrix_Type type;
	Matrix_Shape shape;
	ElementMapping<T> mapping;
	bool touched;
};

template <typename T>
class Matrix_Base: public Matrix_Base_Common<Matrix_Base<T>, T> {
public:
	virtual ~Matrix_Base() {}

	using Matrix_Base_Common<Matrix_Base<T>, T>::copy;
	using Matrix_Base_Common<Matrix_Base<T>, T>::add;

	virtual void copy(const Matrix_Base<T> *in, int rowOffset, int colOffset, int size, int step) =0;
	virtual void add(const T &alpha, const Matrix_Base<T> *a, int rowOffset, int colOffset, int size, int step) =0;
};

template <typename T>
class Matrix_Base<std::complex<T> >: public Matrix_Base_Common<Matrix_Base<std::complex<T> >, std::complex<T> > {
public:
	virtual ~Matrix_Base() {}

	virtual void copyReal(const Matrix_Base<T> *in) =0;
	virtual void copyImag(const Matrix_Base<T> *in) =0;
	virtual void addReal(const T &alpha, const Matrix_Base<T> *a) =0;
	virtual void addImag(const T &alpha, const Matrix_Base<T> *a) =0;
};

}



#endif /* SRC_MATH2_GENERALIZATION_MATRIX_BASE_H_ */
