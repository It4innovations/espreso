
#ifndef SRC_MATH2_GENERALIZATION_MATRIX_DISTRIBUTED_H_
#define SRC_MATH2_GENERALIZATION_MATRIX_DISTRIBUTED_H_

#include "vector_distributed.h"
#include "matrix_base.h"
#include "math2/math2.h"
#include "math2/utils/dofs_distribution.h"

#include <vector>

namespace espreso {

template <template<typename> typename Matrix, typename T>
class Matrix_Distributed: public Matrix_Base<T>
{
public:
	~Matrix_Distributed()
	{

	}

//	Matrix_Distributed* copy()
//	{
//		Matrix_Distributed<Matrix, T> *m = new Matrix_Distributed<Matrix, T>();
//		m->cluster.resize(cluster);
//		return m;
//	}

	Matrix_Distributed* copyPattern()
	{
		Matrix_Distributed<Matrix, T> *m = new Matrix_Distributed<Matrix, T>();
		m->type = m->cluster.type = this->type;
		m->shape = m->cluster.shape = this->shape;
		m->cluster.pattern(cluster);
		return m;
	}

	void store(const char *file)
	{
		math::store(*this, file);
	}

	void fill(const T &value)
	{
		math::fill(cluster, value);
	}

	void fillData(const Matrix_Base<T> *in)
	{
		if (dynamic_cast<const Matrix_Distributed<Matrix, T>*>(in)) {
			math::copy(cluster, static_cast<const Matrix_Distributed<Matrix, T>*>(in)->cluster);
		}
	}

	void scale(const T &alpha)
	{
		math::scale(alpha, cluster);
	}

	void add(const T &alpha, const Matrix_Base<T> *a)
	{
		if (dynamic_cast<const Matrix_Distributed<Matrix, T>*>(a)) {
			math::add(cluster, alpha, static_cast<const Matrix_Distributed<Matrix, T>*>(a)->cluster);
		}
	}

	void add(const T &alpha, const Matrix_Base<T> *a, int rowOffset, int colOffset, int size, int step)
	{
		if (dynamic_cast<const Matrix_Distributed<Matrix, T>*>(a)) {
			math::add(cluster, alpha, static_cast<const Matrix_Distributed<Matrix, T>*>(a)->cluster, rowOffset, colOffset, size, step);
		}
	}

	void sum(const T &alpha, const Matrix_Base<T> *a, const T &beta, const Matrix_Base<T> *b)
	{
		if (dynamic_cast<const Matrix_Distributed<Matrix, T>*>(a) && dynamic_cast<const Matrix_Distributed<Matrix, T>*>(b)) {
			math::sum(cluster, alpha, static_cast<const Matrix_Distributed<Matrix, T>*>(a)->cluster, beta, static_cast<const Matrix_Distributed<Matrix, T>*>(b)->cluster);
		}
	}

	void sum(const T &alpha, const Matrix_Base<T> *a, const T &beta, const Matrix_Base<T> *b,  int rowOffset, int colOffset, int size, int step)
	{
		if (dynamic_cast<const Matrix_Distributed<Matrix, T>*>(a) && dynamic_cast<const Matrix_Distributed<Matrix, T>*>(b)) {
			math::sum(cluster, alpha, static_cast<const Matrix_Distributed<Matrix, T>*>(a)->cluster, beta, static_cast<const Matrix_Distributed<Matrix, T>*>(b)->cluster, rowOffset, colOffset, size, step);
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
};

}

#endif /* SRC_MATH2_GENERALIZATION_MATRIX_DISTRIBUTED_H_ */
