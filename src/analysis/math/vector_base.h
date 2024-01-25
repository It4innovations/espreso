
#ifndef SRC_ANALYSIS_MATH_VECTOR_BASE_H_
#define SRC_ANALYSIS_MATH_VECTOR_BASE_H_

#include "analysis/builder/mapping.h"
#include "math/primitives/vector_dense.h"
#include "math/primitives/vector_sparse.h"
#include <complex>
#include <vector>

namespace espreso {

template <template<typename, typename> typename Vector, typename T> class Vector_Distributed;
template <template<typename, typename> typename Vector, typename T> class Vector_FETI;

template <typename T> class Vector_Base {
public:
	Vector_Base(): filled(false), updated(false) {}
	virtual ~Vector_Base() {};

	virtual void synchronize() =0;

	virtual Vector_Base<T>* copyPattern() =0;
	virtual void store(const char *file) =0;
	virtual void storeTo(std::vector<double> &output) =0;
	virtual void setFrom(std::vector<double> &output) =0;

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
	virtual void copyTo(Vector_FETI<Vector_Dense , T> *a) const =0;
	virtual void copyTo(Vector_FETI<Vector_Sparse, T> *a) const =0;

	virtual void addTo(const T &alpha, Vector_Distributed<Vector_Dense , T> *a) const =0;
	virtual void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, T> *a) const =0;
	virtual void addTo(const T &alpha, Vector_FETI<Vector_Dense , T> *a) const =0;
	virtual void addTo(const T &alpha, Vector_FETI<Vector_Sparse, T> *a) const =0;

	Mapping<T> mapping;
	bool filled, updated;
};

}

#endif /* SRC_ANALYSIS_MATH_VECTOR_BASE_H_ */
