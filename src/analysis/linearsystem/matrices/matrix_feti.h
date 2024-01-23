
#ifndef SRC_MATH2_GENERALIZATION_MATRIX_FETI_H_
#define SRC_MATH2_GENERALIZATION_MATRIX_FETI_H_

#include "vector_feti.h"
#include "matrix_base.h"
#include "math.physics.h"
#include "matrix_feti.decomposition.h"
#include "math/primitives/matrix_dense.h"
#include "math/math.h"

#include <vector>

namespace espreso {

template <template<typename, typename> typename Matrix, typename T>
class Matrix_FETI: public Matrix_Base<T> {
public:
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
		if (spblas.size() == 0) {
			spblas.resize(domains.size());
			#pragma omp parallel for
			for (size_t i = 0; i < spblas.size(); ++i) {
				spblas[i].insert(domains[i]);
			}
		}

		const Vector_FETI<Vector_Dense, T> *_in = dynamic_cast<const Vector_FETI<Vector_Dense, T>*>(in);
		Vector_FETI<Vector_Dense, T> *_out = dynamic_cast<Vector_FETI<Vector_Dense, T>*>(out);
		if (_in && _out) {
			#pragma omp parallel for
			for (size_t i = 0; i < spblas.size(); ++i) {
				spblas[i].apply(_out->domains[i], alpha, beta, _in->domains[i]);
			}
		} else {
			eslog::error("call empty function Matrix_FETI::apply\n");
		}
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

	std::vector<Matrix<T, int> > domains;
	std::vector<SpBLAS<Matrix, T, int> > spblas;
	DOFsDecomposition *decomposition;
};

}

#endif /* SRC_MATH2_GENERALIZATION_MATRIX_FETI_H_ */
