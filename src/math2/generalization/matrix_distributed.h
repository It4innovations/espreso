
#ifndef SRC_MATH2_GENERALIZATION_MATRIX_DISTRIBUTED_H_
#define SRC_MATH2_GENERALIZATION_MATRIX_DISTRIBUTED_H_

#include "vector_distributed.h"
#include "matrix_base.h"
#include "math2/math2.h"
#include "math2/primitives/matrix_dense.h"
#include "math2/generalization/matrix_distributed.h"
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
		in->copyTo(static_cast<Matrix_Distributed<Matrix, T>*>(this));
	}

	void add(const T &alpha, const Matrix_Base<T> *a)
	{
		a->addTo(alpha, static_cast<Matrix_Distributed<Matrix, T>*>(this));
	}

	void apply(const T &alpha, const Vector_Base<T> *in, const T &beta, Vector_Base<T> *out)
	{
		if (dynamic_cast<const Vector_Distributed<Vector_Dense, T>*>(in) && dynamic_cast<const Vector_Distributed<Vector_Dense, T>*>(out)) {
			math::apply<T>(static_cast<Vector_Distributed<Vector_Dense, T>*>(out)->cluster, alpha, cluster, beta, static_cast<const Vector_Distributed<Vector_Dense, T>*>(in)->cluster);
		}
	}

	void copyTo(Matrix_Distributed<Matrix_Dense, T> *a) const
	{
		math::copy(a->cluster, this->cluster);
	}

	void copyTo(Matrix_Distributed<Matrix_CSR, T> *a) const
	{
		math::copy(a->cluster, this->cluster);
	}

	void copyTo(Matrix_Distributed<Matrix_IJV, T> *a) const
	{
		math::copy(a->cluster, this->cluster);
	}

	void addTo(const T &alpha, Matrix_Distributed<Matrix_Dense, T> *a) const
	{
		math::add(a->cluster, alpha, this->cluster);
	}

	void addTo(const T &alpha, Matrix_Distributed<Matrix_CSR, T> *a) const
	{
		math::add(a->cluster, alpha, this->cluster);
	}

	void addTo(const T &alpha, Matrix_Distributed<Matrix_IJV, T> *a) const
	{
		math::add(a->cluster, alpha, this->cluster);
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
		in->copyTo(this, rowOffset, colOffset, size, step);
	}

	void add(const T &alpha, const Matrix_Base<T> *a, int rowOffset, int colOffset, int size, int step)
	{
		a->addTo(alpha, this, rowOffset, colOffset, size, step);
	}

	void copyTo(Matrix_Distributed<Matrix_Dense, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		math::copy(a->cluster, this->cluster, rowOffset, colOffset, size, step);
	}

	void copyTo(Matrix_Distributed<Matrix_CSR, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		math::copy(a->cluster, this->cluster, rowOffset, colOffset, size, step);
	}

	void copyTo(Matrix_Distributed<Matrix_IJV, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		math::copy(a->cluster, this->cluster, rowOffset, colOffset, size, step);
	}

	void addTo(const T &alpha, Matrix_Distributed<Matrix_Dense, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		math::add(a->cluster, alpha, this->cluster, rowOffset, colOffset, size, step);
	}

	void addTo(const T &alpha, Matrix_Distributed<Matrix_CSR, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		math::add(a->cluster, alpha, this->cluster, rowOffset, colOffset, size, step);
	}

	void addTo(const T &alpha, Matrix_Distributed<Matrix_IJV, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		math::add(a->cluster, alpha, this->cluster, rowOffset, colOffset, size, step);
	}

	void copyToReal(Matrix_Distributed<Matrix_Dense, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToReal(Matrix_Distributed<Matrix_CSR, std::complex<T> > *a) const
	{
		math::copy(a->cluster, 0, this->cluster);
	}

	void copyToReal(Matrix_Distributed<Matrix_IJV, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToImag(Matrix_Distributed<Matrix_Dense, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToImag(Matrix_Distributed<Matrix_CSR, std::complex<T> > *a) const
	{
		math::copy(a->cluster, 1, this->cluster);
	}

	void copyToImag(Matrix_Distributed<Matrix_IJV, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void addToReal(const T &alpha, Matrix_Distributed<Matrix_Dense, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void addToReal(const T &alpha, Matrix_Distributed<Matrix_CSR  , std::complex<T> > *a) const
	{
		math::add(a->cluster, 0, alpha, this->cluster);
	}

	void addToReal(const T &alpha, Matrix_Distributed<Matrix_IJV  , std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void addToImag(const T &alpha, Matrix_Distributed<Matrix_Dense, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void addToImag(const T &alpha, Matrix_Distributed<Matrix_CSR  , std::complex<T> > *a) const
	{
		math::add(a->cluster, 1, alpha, this->cluster);
	}

	void addToImag(const T &alpha, Matrix_Distributed<Matrix_IJV  , std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}
};

template <template<typename> typename Matrix, typename T>
class Matrix_Distributed<Matrix, std::complex<T> >: public Matrix_Distributed_Common<Matrix, std::complex<T> > {
public:
	void copyReal(const Matrix_Base<T> *in)
	{
		in->copyToReal(this);
	}

	void copyImag(const Matrix_Base<T> *in)
	{
		in->copyToImag(this);
	}

	void addReal(const T &alpha, const Matrix_Base<T> *a)
	{
		a->addToReal(alpha, this);
	}

	void addImag(const T &alpha, const Matrix_Base<T> *a)
	{
		a->addToImag(alpha, this);
	}
};

}

#endif /* SRC_MATH2_GENERALIZATION_MATRIX_DISTRIBUTED_H_ */
