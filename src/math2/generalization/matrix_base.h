
#ifndef SRC_MATH2_GENERALIZATION_MATRIX_BASE_H_
#define SRC_MATH2_GENERALIZATION_MATRIX_BASE_H_

#include "vector_base.h"

#include "analysis/composer/elementmapping.h"

#include <complex>

namespace espreso {

template <typename T> class Matrix_Base;

template <typename T>
class Matrix_Base_Common {
public:
	Matrix_Base_Common(Matrix_Type type): type(type), shape(Matrix_Shape::FULL), touched(false) {}

	virtual ~Matrix_Base_Common() {};

	virtual void commit() =0;

	virtual Matrix_Base<T>* copyPattern() =0;
	virtual void store(const char *file) =0;

	virtual void set(const T &value) =0;
	virtual void scale(const T &alpha) =0;

	virtual void copy(const Matrix_Base<T> *in) =0;
	virtual void add(const T &alpha, const Matrix_Base<T> *a) =0;
	virtual void apply(const T &alpha, const Vector_Base<T> *in, const T &beta, Vector_Base<T> *out) =0;

	Matrix_Type type;
	Matrix_Shape shape;
	ElementMapping<T> mapping;
	bool touched;
};

template <typename T>
class Matrix_Base: public Matrix_Base_Common<T> {
public:
	Matrix_Base(): Matrix_Base_Common<T>(Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC) {}
	virtual ~Matrix_Base() {}

	using Matrix_Base_Common<T>::copy;
	using Matrix_Base_Common<T>::add;

	virtual void copy(const Matrix_Base<T> *in, int rowOffset, int colOffset, int size, int step) =0;
	virtual void add(const T &alpha, const Matrix_Base<T> *a, int rowOffset, int colOffset, int size, int step) =0;
};

template <typename T>
class Matrix_Base<std::complex<T> >: public Matrix_Base_Common<std::complex<T> > {
public:
	Matrix_Base(): Matrix_Base_Common<std::complex<T> >(Matrix_Type::COMPLEX_STRUCTURALLY_SYMMETRIC) {}
	virtual ~Matrix_Base() {}

	virtual void copyReal(const Matrix_Base<T> *in) =0;
	virtual void copyImag(const Matrix_Base<T> *in) =0;
	virtual void addReal(const T &alpha, const Matrix_Base<T> *a) =0;
	virtual void addImag(const T &alpha, const Matrix_Base<T> *a) =0;
};

}



#endif /* SRC_MATH2_GENERALIZATION_MATRIX_BASE_H_ */
