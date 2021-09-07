
#ifndef SRC_MATH2_GENERALIZATION_MATRIX_BASE_H_
#define SRC_MATH2_GENERALIZATION_MATRIX_BASE_H_

#include "vector_base.h"

#include "analysis/composer/elementmapping.h"

namespace espreso {

template <typename T>
class Matrix_Base
{
public:
	Matrix_Base(): type(Matrix_Type::REAL_UNSYMMETRIC), shape(Matrix_Shape::FULL), touched(false) {}

	virtual ~Matrix_Base() {};

	virtual void commit() =0;

//	virtual Matrix_Base* copy() =0;
	virtual Matrix_Base* copyPattern() =0;
	virtual void store(const char *file) =0;

	virtual void fill(const T &value) =0;
	virtual void fillData(const Matrix_Base *in) =0;

	virtual void scale(const T &alpha) =0;
	virtual void add(const T &alpha, const Matrix_Base *a) =0;
	virtual void add(const T &alpha, const Matrix_Base *a, int rowOffset, int colOffset, int size, int step) =0;
	virtual void sum(const T &alpha, const Matrix_Base *a, const T &beta, const Matrix_Base *b) =0;
	virtual void sum(const T &alpha, const Matrix_Base *a, const T &beta, const Matrix_Base *b, int rowOffset, int colOffset, int size, int step) =0;

	virtual void apply(const T &alpha, const Vector_Base<T> *in, const T &beta, Vector_Base<T> *out) =0;

	Matrix_Type type;
	Matrix_Shape shape;
	ElementMapping<T> mapping;
	bool touched;
};

}



#endif /* SRC_MATH2_GENERALIZATION_MATRIX_BASE_H_ */
