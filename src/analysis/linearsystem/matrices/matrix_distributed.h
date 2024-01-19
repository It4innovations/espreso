
#ifndef SRC_MATH2_GENERALIZATION_MATRIX_DISTRIBUTED_H_
#define SRC_MATH2_GENERALIZATION_MATRIX_DISTRIBUTED_H_

#include "vector_distributed.h"
#include "matrix_base.h"
#include "matrix_distributed.distribution.h"
#include "math.physics.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "matrix_distributed.apply.h"
#include "matrix_distributed.synchronization.h"

#include <vector>

namespace espreso {

template <template<typename, typename, template<typename> typename> typename Matrix, typename T>
class Matrix_Distributed: public Matrix_Base<T> {
public:
	void commit()
	{
		applyData.commit(*static_cast<Matrix_Distributed<Matrix, T>*>(this));
	}

	void initApply()
	{
		applyData.init(*static_cast<Matrix_Distributed<Matrix, T>*>(this));
	}

	void synchronize()
	{
		synchronization->gatherFromUpper(*static_cast<Matrix_Distributed<Matrix, T>*>(this));
	}

	Matrix_Base<T>* copyPattern()
	{
		Matrix_Distributed<Matrix, T> *m = new Matrix_Distributed<Matrix, T>();
		m->type = m->cluster.type = this->type;
		m->shape = m->cluster.shape = this->shape;
		m->cluster.pattern(cluster);
		m->distribution = this->distribution;
		m->applyData = this->applyData;
		m->synchronization = this->synchronization;
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
		this->touched = true;
		in->copyTo(static_cast<Matrix_Distributed<Matrix, T>*>(this));
	}

	void add(const T &alpha, const Matrix_Base<T> *a)
	{
		a->addTo(alpha, static_cast<Matrix_Distributed<Matrix, T>*>(this));
	}

	void apply(const T &alpha, const Vector_Base<T> *in, const T &beta, Vector_Base<T> *out)
	{
		if (dynamic_cast<const Vector_Distributed<Vector_Dense, T>*>(in) && dynamic_cast<const Vector_Distributed<Vector_Dense, T>*>(out)) {
			applyData.apply(static_cast<Vector_Distributed<Vector_Dense, T>*>(out), alpha, beta, static_cast<const Vector_Distributed<Vector_Dense, T>*>(in));
		}
	}

	void copyTo(Matrix_Distributed<Matrix_CSR, T> *a) const
	{
		math::copy(a->cluster, this->cluster);
	}

	void copyTo(Matrix_FETI<Matrix_CSR, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void addTo(const T &alpha, Matrix_Distributed<Matrix_CSR, T> *a) const
	{
		math::add(a->cluster, alpha, this->cluster);
	}

	void addTo(const T &alpha, Matrix_FETI<Matrix_CSR, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	Matrix<T, esint, cpu_allocator> cluster;
	DOFsDistribution *distribution;
	Data_Apply<Matrix, T> applyData;
	Data_Synchronization<Matrix, T> *synchronization;
};

}

#endif /* SRC_MATH2_GENERALIZATION_MATRIX_DISTRIBUTED_H_ */
