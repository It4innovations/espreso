
#ifndef SRC_MATH2_GENERALIZATION_MATRIX_FETI_H_
#define SRC_MATH2_GENERALIZATION_MATRIX_FETI_H_

#include "vector_feti.h"
#include "matrix_base.h"
#include "math2/math2.h"

#include <vector>

namespace espreso {

template <template<typename> typename Matrix, typename T>
class Matrix_FETI: public Matrix_Base<T>
{
public:
	~Matrix_FETI()
	{

	}

//	Matrix_FETI* copy()
//	{
//		Matrix_FETI<Matrix, T> *m = new Matrix_FETI<Matrix, T>();
//		m->domains.resize(domains.size());
//		for (size_t d = 0; d < domains.size(); ++d) {
//			m->domains[d].resize(domains[d]);
//		}
//		return m;
//	}

	void commit()
	{
		for (size_t d = 0; d < domains.size(); ++d) {
			math::commit(domains[d]);
		}
	}

	Matrix_FETI* copyPattern()
	{
		Matrix_FETI<Matrix, T> *m = new Matrix_FETI<Matrix, T>();
		m->domains.resize(domains.size());
		for (size_t d = 0; d < domains.size(); ++d) {
			m->type = m->domains[d].type = this->type;
			m->shape = m->domains[d].shape = this->shape;
			m->domains[d].pattern(domains[d]);
		}
		return m;
	}

	void store(const char *file)
	{
		math::store(*this, file);
	}

	void fill(const T &value)
	{
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); ++d) {
			math::fill(domains[d], value);
		}
	}

	void fillData(const Matrix_Base<T> *in)
	{
		if (dynamic_cast<const Matrix_FETI<Matrix, T>*>(in)) {
			#pragma omp parallel for
			for (size_t d = 0; d < domains.size(); ++d) {
				math::copy(domains[d], static_cast<const Matrix_FETI<Matrix, T>*>(in)->domains[d]);
			}
		}
	}

	void scale(const T &alpha)
	{
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); ++d) {
			math::scale(alpha, domains[d]);
		}
	}

	void add(const T &alpha, const Matrix_Base<T> *a)
	{
		if (dynamic_cast<const Matrix_FETI<Matrix, T>*>(a)) {
			#pragma omp parallel for
			for (size_t d = 0; d < domains.size(); ++d) {
				math::add(domains[d], alpha, static_cast<const Matrix_FETI<Matrix, T>*>(a)->domains[d]);
			}
		}
	}

	void add(const T &alpha, const Matrix_Base<T> *a, int rowOffset, int colOffset, int size, int step)
	{
		if (dynamic_cast<const Matrix_FETI<Matrix, T>*>(a)) {
			#pragma omp parallel for
			for (size_t d = 0; d < domains.size(); ++d) {
				math::add(domains[d], alpha, static_cast<const Matrix_FETI<Matrix, T>*>(a)->domains[d], rowOffset, colOffset, size, step);
			}
		}
	}

	void sum(const double &alpha, const Matrix_Base<double> *a, const double &beta, const Matrix_Base<double> *b)
	{
		if (dynamic_cast<const Matrix_FETI<Matrix, T>*>(a) && dynamic_cast<const Matrix_FETI<Matrix, T>*>(b)) {
			#pragma omp parallel for
			for (size_t d = 0; d < domains.size(); ++d) {
				math::sum(domains[d], alpha, static_cast<const Matrix_FETI<Matrix, T>*>(a)->domains[d], beta, static_cast<const Matrix_FETI<Matrix, T>*>(b)->domains[d]);
			}
		}
	}

	void sum(const std::complex<double> &alpha, const Matrix_Base<std::complex<double>> *a, const std::complex<double> &beta, const Matrix_Base<std::complex<double>> *b)
	{
		if (dynamic_cast<const Matrix_FETI<Matrix, T>*>(a) && dynamic_cast<const Matrix_FETI<Matrix, T>*>(b)) {
			#pragma omp parallel for
			for (size_t d = 0; d < domains.size(); ++d) {
				math::sum(domains[d], alpha, static_cast<const Matrix_FETI<Matrix, T>*>(a)->domains[d], beta, static_cast<const Matrix_FETI<Matrix, T>*>(b)->domains[d]);
			}
		}
	}


	void sum(const T &alpha, const Matrix_Base<T> *a, const T &beta, const Matrix_Base<T> *b, int rowOffset, int colOffset, int size, int step)
	{
		if (dynamic_cast<const Matrix_FETI<Matrix, T>*>(a) && dynamic_cast<const Matrix_FETI<Matrix, T>*>(b)) {
			#pragma omp parallel for
			for (size_t d = 0; d < domains.size(); ++d) {
				math::sum(domains[d], alpha, static_cast<const Matrix_FETI<Matrix, T>*>(a)->domains[d], beta, static_cast<const Matrix_FETI<Matrix, T>*>(b)->domains[d], rowOffset, colOffset, size, step);
			}
		}
	}

	void apply(const T &alpha, const Vector_Base<T> *in, const T &beta, Vector_Base<T> *out)
	{
		if (dynamic_cast<const Vector_FETI<Vector_Dense, T>*>(in) && dynamic_cast<const Vector_FETI<Vector_Dense, T>*>(out)) {
			#pragma omp parallel for
			for (size_t d = 0; d < domains.size(); ++d) {
				math::apply<T>(static_cast<Vector_FETI<Vector_Dense, T>*>(out)->domains[d], alpha, domains[d], beta, static_cast<const Vector_FETI<Vector_Dense, T>*>(in)->domains[d]);
			}
		}
	}

	std::vector<Matrix<T> > domains;
};

}

#endif /* SRC_MATH2_GENERALIZATION_MATRIX_FETI_H_ */
