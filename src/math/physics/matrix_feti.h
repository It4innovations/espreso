
#ifndef SRC_MATH2_GENERALIZATION_MATRIX_FETI_H_
#define SRC_MATH2_GENERALIZATION_MATRIX_FETI_H_


#include "vector_feti.h"
#include "matrix_base.h"
#include "math/physics/math.physics.h"
#include "math/primitives/matrix_dense.h"
#include "math/utils/feti/decomposition.h"

#include <vector>

namespace espreso {

template <template<typename, typename, template<typename> typename> typename Matrix, typename T>
class Matrix_FETI: public Matrix_Base<T> {
public:
	void commit()
	{
//		spblas.resize(this->domains.size());
//		#pragma omp parallel for
//		for (size_t d = 0; d < this->domains.size(); ++d) {
//			spblas.commit(this->domains[d]);
//		}
	}

	void combine(const Matrix_FETI<Matrix, T> &A, const Matrix_FETI<Matrix, T> &B)
	{
		for (size_t d = 0; d < this->domains.size(); ++d) {
			if (A.domains[d].type != B.domains[d].type) {
				eslog::error("cannot combine matrices of different types.\n");
			}
			if (A.domains[d].shape != B.domains[d].shape) {
				eslog::error("cannot combine matrices of different shapes.\n");
			}
		}
		this->type = A.type;
		this->shape = A.shape;
		this->domains.resize(A.domains.size());
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			this->domains[d].type = A.domains[d].type;
			this->domains[d].shape = A.domains[d].shape;
			math::combine(this->domains[d], A.domains[d], B.domains[d]);
		}
	}

	void sumCombined(const T &alpha, const Matrix_FETI<Matrix, T> &A, const Matrix_FETI<Matrix, T> &B)
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::sumCombined(this->domains[d], alpha, A.domains[d], B.domains[d]);
		}
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
		this->touched = true;
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

	void copyTo(Matrix_Distributed<Matrix_CSR, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyTo(Matrix_FETI<Matrix_CSR, T> *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], this->domains[d]);
		}
	}

	void addTo(const T &alpha, Matrix_Distributed<Matrix_CSR, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void addTo(const T &alpha, Matrix_FETI<Matrix_CSR, T> *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::add(a->domains[d], alpha, this->domains[d]);
		}
	}

	std::vector<Matrix<T, int, cpu_allocator> > domains;
	DOFsDecomposition *decomposition;
};

}

#endif /* SRC_MATH2_GENERALIZATION_MATRIX_FETI_H_ */
