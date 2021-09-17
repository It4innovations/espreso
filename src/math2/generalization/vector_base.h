
#ifndef SRC_MATH2_GENERALIZATION_VECTOR_BASE_H_
#define SRC_MATH2_GENERALIZATION_VECTOR_BASE_H_

#include "analysis/composer/elementmapping.h"
#include "math2/primitives/vector_sparse.h"

#include <complex>
#include <vector>

namespace espreso {

template <typename T> class Vector_Dense;
template <template<typename> typename Vector, typename T> class Vector_Distributed;

template <typename T>
class Vector_Base
{
public:
	Vector_Base(): touched(false) {}
	virtual ~Vector_Base() {};

//	virtual Vector_Base* copy() =0;
	virtual Vector_Base* copyPattern() =0;
	virtual void store(const char *file) =0;
	virtual void store(std::vector<T> &output) =0;

	virtual void fill(const T &value) =0;
	virtual void fillData(const Vector_Base *in) =0;
	virtual void fillData(const Vector_Base *in, int offset, int size, int step) =0;

	virtual void scale(const T &alpha) =0;
	//virtual void add(const T &alpha, const Vector_Base *a) =0;
	//virtual void add(const double &alpha, const Vector_Base<std::complex<double>> *a) =0;
	virtual void add(const double &alpha, const Vector_Base<double> *a) =0;
	virtual void add(const T &alpha, const Vector_Base *a, int offset, int size, int step) =0;
	virtual void add_imag(const double &alpha, const Vector_Base<double> *a) =0;

	virtual void sum(const T &alpha, const Vector_Base *a, const T &beta, const Vector_Base *b) =0;
	virtual void sum(const T &alpha, const Vector_Base *a, const T &beta, const Vector_Base *b, int offset, int size, int step) =0;

	virtual T norm() =0;
	virtual T max() =0;
	virtual T absmax() =0;
	virtual T dot(const Vector_Base *other) =0;

	ElementMapping<T> mapping;
	bool touched;

	virtual void addTo(const T &alpha, Vector_Distributed<Vector_Dense, T> *a) const =0;
	virtual void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, T> *a) const =0;
	virtual void addTo(const T &alpha, Vector_Distributed<Vector_Dense, T> *a, int offset, int size, int step) const =0;
	virtual void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, T> *a, int offset, int size, int step) const =0;
};

}

#endif /* SRC_MATH2_GENERALIZATION_VECTOR_BASE_H_ */
