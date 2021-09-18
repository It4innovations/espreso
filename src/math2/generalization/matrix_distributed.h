
#ifndef SRC_MATH2_GENERALIZATION_MATRIX_DISTRIBUTED_H_
#define SRC_MATH2_GENERALIZATION_MATRIX_DISTRIBUTED_H_

#include "vector_distributed.h"
#include "matrix_base.h"
#include "math2/math2.h"
#include "math2/utils/dofs_distribution.h"
#include "math2/utils/utils_distributed.h"

#include <vector>

namespace espreso {

template <template<typename> typename Matrix, typename T> class Matrix_Distributed;

template <template<typename> typename Matrix, typename T>
class Matrix_Distributed_Common: public Matrix_Base<T> {
public:
	void commit()
	{
		math::commit(cluster);
	}

	Matrix_Base<T>* copyPattern()
	{
		Matrix_Distributed<Matrix, T> *m = new Matrix_Distributed<Matrix, T>();
		m->type = m->cluster.type = this->type;
		m->shape = m->cluster.shape = this->shape;
		m->cluster.pattern(cluster);
		return m;
	}

	void store(const char *file)
	{
		math::store(*static_cast<Matrix_Distributed<Matrix, T>*>(this), file);
	}

	void set(const T &value)
	{
		math::set(cluster, value);
	}

	void scale(const T &alpha)
	{
		math::scale(alpha, cluster);
	}

	void copy(const Matrix_Base<T> *in)
	{
		if (dynamic_cast<const Matrix_Distributed<Matrix, T>*>(in)) {
			math::copy(cluster, static_cast<const Matrix_Distributed<Matrix, T>*>(in)->cluster);
		}
	}

	void add(const T &alpha, const Matrix_Base<T> *a)
	{
		if (dynamic_cast<const Matrix_Distributed<Matrix, T>*>(a)) {
			math::add(cluster, alpha, static_cast<const Matrix_Distributed<Matrix, T>*>(a)->cluster);
		}
	}

	void apply(const T &alpha, const Vector_Base<T> *in, const T &beta, Vector_Base<T> *out)
	{
		if (dynamic_cast<const Vector_Distributed<Vector_Dense, T>*>(in) && dynamic_cast<const Vector_Distributed<Vector_Dense, T>*>(out)) {
			math::apply<T>(static_cast<Vector_Distributed<Vector_Dense, T>*>(out)->cluster, alpha, cluster, beta, static_cast<const Vector_Distributed<Vector_Dense, T>*>(in)->cluster);
		}
	}

	Matrix<T> cluster;
	DOFsDistribution distribution;
	DataSynchronization synchronization;

};

template <template<typename> typename Matrix, typename T>
class Matrix_Distributed: public Matrix_Distributed_Common<Matrix, T> {
public:
	void copy(const Matrix_Base<T> *in, int rowOffset, int colOffset, int size, int step)
	{
		if (dynamic_cast<const Matrix_Distributed<Matrix, T>*>(in)) {
			math::copy(this->cluster, static_cast<const Matrix_Distributed<Matrix, T>*>(in)->cluster, rowOffset, colOffset, size, step);
		}
	}

	void add(const T &alpha, const Matrix_Base<T> *a, int rowOffset, int colOffset, int size, int step)
	{
		if (dynamic_cast<const Matrix_Distributed<Matrix, T>*>(a)) {
			math::add(this->cluster, alpha, static_cast<const Matrix_Distributed<Matrix, T>*>(a)->cluster, rowOffset, colOffset, size, step);
		}
	}
};

template <template<typename> typename Matrix, typename T>
class Matrix_Distributed<Matrix, std::complex<T> >: public Matrix_Distributed_Common<Matrix, std::complex<T> > {
public:
	void copyReal(const Matrix_Base<T> *in)
	{
		if (dynamic_cast<const Matrix_Distributed<Matrix, T>*>(in)) {
			math::copy(this->cluster, 0, static_cast<const Matrix_Distributed<Matrix, T>*>(in)->cluster);
		}
	}

	void copyImag(const Matrix_Base<T> *in)
	{
		if (dynamic_cast<const Matrix_Distributed<Matrix, T>*>(in)) {
			math::copy(this->cluster, 1, static_cast<const Matrix_Distributed<Matrix, T>*>(in)->cluster);
		}
	}

	void addReal(const T &alpha, const Matrix_Base<T> *a)
	{
		if (dynamic_cast<const Matrix_Distributed<Matrix, T>*>(a)) {
			math::add(this->cluster, 0, alpha, static_cast<const Matrix_Distributed<Matrix, T>*>(a)->cluster);
		}
	}

	void addImag(const T &alpha, const Matrix_Base<T> *a)
	{
		if (dynamic_cast<const Matrix_Distributed<Matrix, T>*>(a)) {
			math::add(this->cluster, 1, alpha, static_cast<const Matrix_Distributed<Matrix, T>*>(a)->cluster);
		}
	}
};

}

#endif /* SRC_MATH2_GENERALIZATION_MATRIX_DISTRIBUTED_H_ */
