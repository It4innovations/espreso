
#ifndef SRC_MATH2_GENERALIZATION_MATRIX_FETI_H_
#define SRC_MATH2_GENERALIZATION_MATRIX_FETI_H_

#include "vector_feti.h"
#include "matrix_base.h"
#include "math2/math2.h"
#include "math2/primitives/matrix_dense.h"
#include "math2/utils/decomposed/decomposition.h"

#include <vector>

namespace espreso {

template <template<typename> typename Matrix, typename T> class Matrix_FETI;

template <template<typename> typename Matrix, typename T>
class Matrix_FETI_Common: public Matrix_Base<T> {
public:
	void commit()
	{
//		applyData.commit(*static_cast<Matrix_FETI<Matrix, T>*>(this));
	}

	void initApply()
	{
//		applyData.init(*static_cast<Matrix_FETI<Matrix, T>*>(this));
	}

	void synchronize()
	{
//		synchronization->gatherFromUpper(*static_cast<Matrix_FETI<Matrix, T>*>(this));
	}

	Matrix_Base<T>* copyPattern()
	{
		Matrix_FETI<Matrix, T> *m = new Matrix_FETI<Matrix, T>();
		m->type = this->type;
		m->shape = this->shape;
		m->domains.resize(domains.size());
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); ++d) {
			m->domains[d].type = domains[d].type;
			m->domains[d].shape = domains[d].shape;
			m->domains[d].pattern(domains[d]);
		}
		return m;
	}

	void store(const char *file)
	{
		math::store(*static_cast<Matrix_FETI<Matrix, T>*>(this), file);
	}

	void set(const T &value)
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::set(this->domains[d], value);
		}
	}

	void scale(const T &alpha)
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::scale(alpha, this->domains[d]);
		}
	}

	void copy(const Matrix_Base<T> *in)
	{
		in->copyTo(static_cast<Matrix_FETI<Matrix, T>*>(this));
	}

	void add(const T &alpha, const Matrix_Base<T> *a)
	{
		a->addTo(alpha, static_cast<Matrix_FETI<Matrix, T>*>(this));
	}

	void apply(const T &alpha, const Vector_Base<T> *in, const T &beta, Vector_Base<T> *out)
	{
		eslog::error("call empty function\n");
	}

	void copyTo(Matrix_Distributed<Matrix_Dense, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyTo(Matrix_Distributed<Matrix_CSR, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyTo(Matrix_Distributed<Matrix_IJV, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyTo(Matrix_FETI<Matrix_Dense, T> *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], this->domains[d]);
		}
	}

	void copyTo(Matrix_FETI<Matrix_CSR, T> *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], this->domains[d]);
		}
	}

	void copyTo(Matrix_FETI<Matrix_IJV, T> *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], this->domains[d]);
		}
	}

	void addTo(const T &alpha, Matrix_Distributed<Matrix_Dense, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void addTo(const T &alpha, Matrix_Distributed<Matrix_CSR, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void addTo(const T &alpha, Matrix_Distributed<Matrix_IJV, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void addTo(const T &alpha, Matrix_FETI<Matrix_Dense, T> *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::add(a->domains[d], alpha, this->domains[d]);
		}
	}

	void addTo(const T &alpha, Matrix_FETI<Matrix_CSR, T> *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::add(a->domains[d], alpha, this->domains[d]);
		}
	}

	void addTo(const T &alpha, Matrix_FETI<Matrix_IJV, T> *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::add(a->domains[d], alpha, this->domains[d]);
		}
	}

	std::vector<Matrix<T> > domains;
	DOFsDecomposition *decomposition;
};

template <template<typename> typename Matrix, typename T>
class Matrix_FETI: public Matrix_FETI_Common<Matrix, T> {
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
		eslog::error("call empty function\n");
	}

	void copyTo(Matrix_Distributed<Matrix_CSR, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		eslog::error("call empty function\n");
	}

	void copyTo(Matrix_Distributed<Matrix_IJV, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		eslog::error("call empty function\n");
	}

	void copyTo(Matrix_FETI<Matrix_Dense, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], this->domains[d], rowOffset, colOffset, size, step);
		}
	}

	void copyTo(Matrix_FETI<Matrix_CSR, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], this->domains[d], rowOffset, colOffset, size, step);
		}
	}

	void copyTo(Matrix_FETI<Matrix_IJV, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], this->domains[d], rowOffset, colOffset, size, step);
		}
	}

	void addTo(const T &alpha, Matrix_Distributed<Matrix_Dense, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		eslog::error("call empty function\n");
	}

	void addTo(const T &alpha, Matrix_Distributed<Matrix_CSR, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		eslog::error("call empty function\n");
	}

	void addTo(const T &alpha, Matrix_Distributed<Matrix_IJV, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		eslog::error("call empty function\n");
	}

	void addTo(const T &alpha, Matrix_FETI<Matrix_Dense, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::add(a->domains[d], alpha, this->domains[d], rowOffset, colOffset, size, step);
		}
	}

	void addTo(const T &alpha, Matrix_FETI<Matrix_CSR, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::add(a->domains[d], alpha, this->domains[d], rowOffset, colOffset, size, step);
		}
	}

	void addTo(const T &alpha, Matrix_FETI<Matrix_IJV, T> *a, int rowOffset, int colOffset, int size, int step) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::add(a->domains[d], alpha, this->domains[d], rowOffset, colOffset, size, step);
		}
	}

	void copyToReal(Matrix_Distributed<Matrix_Dense, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToReal(Matrix_Distributed<Matrix_CSR, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToReal(Matrix_Distributed<Matrix_IJV, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToReal(Matrix_FETI<Matrix_Dense, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToReal(Matrix_FETI<Matrix_CSR, std::complex<T> > *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], 0, this->domains[d]);
		}
	}

	void copyToReal(Matrix_FETI<Matrix_IJV, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToImag(Matrix_Distributed<Matrix_Dense, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToImag(Matrix_Distributed<Matrix_CSR, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToImag(Matrix_Distributed<Matrix_IJV, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToImag(Matrix_FETI<Matrix_Dense, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToImag(Matrix_FETI<Matrix_CSR, std::complex<T> > *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], 1, this->domains[d]);
		}
	}

	void copyToImag(Matrix_FETI<Matrix_IJV, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void addToReal(const T &alpha, Matrix_Distributed<Matrix_Dense, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void addToReal(const T &alpha, Matrix_Distributed<Matrix_CSR  , std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void addToReal(const T &alpha, Matrix_Distributed<Matrix_IJV  , std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void addToReal(const T &alpha, Matrix_FETI<Matrix_Dense, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void addToReal(const T &alpha, Matrix_FETI<Matrix_CSR  , std::complex<T> > *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::add(a->domains[d], 0, alpha, this->domains[d]);
		}
	}

	void addToReal(const T &alpha, Matrix_FETI<Matrix_IJV  , std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void addToImag(const T &alpha, Matrix_Distributed<Matrix_Dense, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void addToImag(const T &alpha, Matrix_Distributed<Matrix_CSR  , std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void addToImag(const T &alpha, Matrix_Distributed<Matrix_IJV  , std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void addToImag(const T &alpha, Matrix_FETI<Matrix_Dense, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void addToImag(const T &alpha, Matrix_FETI<Matrix_CSR  , std::complex<T> > *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::add(a->domains[d], 1, alpha, this->domains[d]);
		}
	}

	void addToImag(const T &alpha, Matrix_FETI<Matrix_IJV  , std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}
};

template <template<typename> typename Matrix, typename T>
class Matrix_FETI<Matrix, std::complex<T> >: public Matrix_FETI_Common<Matrix, std::complex<T> > {
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

#endif /* SRC_MATH2_GENERALIZATION_MATRIX_FETI_H_ */
