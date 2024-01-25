
#ifndef SRC_ANALYSIS_MATH_MATRIX_DISTRIBUTED_H_
#define SRC_ANALYSIS_MATH_MATRIX_DISTRIBUTED_H_

#include "analysis/math/math.physics.h"
#include "analysis/math/matrix_base.h"
#include "analysis/math/vector_distributed.h"
#include "analysis/builder/direct.apply.h"
#include "analysis/builder/direct.synchronization.h"
#include "analysis/builder/direct.decomposition.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"

#include <vector>

namespace espreso {

template <typename T>
class Matrix_Distributed: public Matrix_Base<T> {
public:
	void synchronize()
	{
		_sync->gatherFromUpper(*static_cast<Matrix_Distributed<T>*>(this));
	}

	Matrix_Base<T>* copyPattern()
	{
		Matrix_Distributed<T> *m = new Matrix_Distributed<T>();
		m->type = m->cluster.type = this->type;
		m->shape = m->cluster.shape = this->shape;
		m->cluster.pattern(cluster);
		m->decomposition = this->decomposition;
		m->_apply = this->_apply;
		m->_sync = this->_sync;
		return m;
	}

	void store(const char *file)
	{
		math::store(*static_cast<Matrix_Distributed<T>*>(this), file);
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
		in->copyTo(static_cast<Matrix_Distributed<T>*>(this));
	}

	void add(const T &alpha, const Matrix_Base<T> *a)
	{
		a->addTo(alpha, static_cast<Matrix_Distributed<T>*>(this));
	}

	void apply(const T &alpha, const Vector_Base<T> *in, const T &beta, Vector_Base<T> *out)
	{
		if (dynamic_cast<const Vector_Distributed<Vector_Dense, T>*>(in) && dynamic_cast<const Vector_Distributed<Vector_Dense, T>*>(out)) {
			_apply->apply(*this, static_cast<Vector_Distributed<Vector_Dense, T>*>(out), alpha, beta, static_cast<const Vector_Distributed<Vector_Dense, T>*>(in));
		}
	}

	void copyTo(Matrix_Distributed<T> *a) const
	{
		math::copy(a->cluster, this->cluster);
	}

	void copyTo(Matrix_FETI<T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void addTo(const T &alpha, Matrix_Distributed<T> *a) const
	{
		math::add(a->cluster, alpha, this->cluster);
	}

	void addTo(const T &alpha, Matrix_FETI<T> *a) const
	{
		eslog::error("call empty function\n");
	}

	Matrix_CSR<T, esint> cluster;
	DirectDecomposition *decomposition;

	Matrix_CSR_Apply<T> *_apply;
	Matrix_CSR_Sync<T> *_sync;
};

}

#endif /* SRC_ANALYSIS_MATH_MATRIX_DISTRIBUTED_H_ */
