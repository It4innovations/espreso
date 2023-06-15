
#ifndef SRC_MATH2_GENERALIZATION_VECTOR_FETI_H_
#define SRC_MATH2_GENERALIZATION_VECTOR_FETI_H_

#include <math/utils/feti/decomposition.h>
#include "vector_base.h"
#include "esinfo/eslog.h"
#include "math/physics/math.physics.h"
#include <vector>

namespace espreso {

template <template<typename> typename Vector, typename T>
class Vector_FETI_Common: public Vector_Base<T> {
public:
	void synchronize()
	{
//		synchronization->synchronize(*static_cast<Vector_FETI<Vector, T>*>(this));
	}

	void scatter()
	{
//		synchronization->scatterToUpper(*static_cast<Vector_FETI<Vector, T>*>(this));
	}

	Vector_Base<T>* copyPattern()
	{
		Vector_FETI<Vector, T> *m = new Vector_FETI<Vector, T>();
		m->domains.resize(domains.size());
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); ++d) {
			m->domains[d].pattern(domains[d]);
		}
		m->decomposition = decomposition;
		return m;
	}

	void store(const char *file)
	{
		math::store(*static_cast<Vector_FETI<Vector, T>*>(this), file);
	}

	void set(const T &value)
	{
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); ++d) {
			math::set(domains[d], value);
		}
	}

	void scale(const T &alpha)
	{
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); ++d) {
			math::scale(alpha, domains[d]);
		}
	}

	void copy(const Vector_Base<T> *in)
	{
		this->touched = true;
		in->copyTo(static_cast<Vector_FETI<Vector, T>*>(this));
	}

	void add(const T &alpha, const Vector_Base<T> *a)
	{
		a->addTo(alpha, static_cast<Vector_FETI<Vector, T>*>(this));
	}

	T norm()
	{
//		T dot = math::dot(cluster.size - distribution->halo.size(), cluster.vals + distribution->halo.size(), 1, cluster.vals + distribution->halo.size(), 1);
//		Communication::allReduce(&dot, NULL, 1, MPI_DOUBLE, MPI_SUM);
//		return std::sqrt(dot);
		eslog::error("call empty function: max\n");
		return T{};
	}

	T max()
	{
		eslog::error("call empty function: max\n");
		return T{};
	}

	T absmax()
	{
		eslog::error("call empty function: absmax\n");
		return T{};
	}

	T dot(const Vector_Base<T> *other)
	{
		eslog::error("call empty function: dot\n");
		return T{};
	}

	void copyTo(Vector_Distributed<Vector_Dense , T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyTo(Vector_Distributed<Vector_Sparse, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyTo(Vector_FETI<Vector_Dense , T> *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); ++d) {
			math::copy(a->domains[d], domains[d]);
		}
	}

	void copyTo(Vector_FETI<Vector_Sparse, T> *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); ++d) {
			math::copy(a->domains[d], domains[d]);
		}
	}

	void addTo(const T &alpha, Vector_Distributed<Vector_Dense, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void addTo(const T &alpha, Vector_FETI<Vector_Dense, T> *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); ++d) {
			math::add(a->domains[d], alpha, domains[d]);
		}
	}

	void addTo(const T &alpha, Vector_FETI<Vector_Sparse, T> *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); ++d) {
			math::add(a->domains[d], alpha, domains[d]);
		}
	}

	std::vector<Vector<T> > domains;
	DOFsDecomposition *decomposition;
};

template <template<typename> typename Vector, typename T>
class Vector_FETI: public Vector_FETI_Common<Vector, T> {

	void _store(const std::vector<Vector_Dense<double> > &domains, const DOFsDecomposition *decomposition, std::vector<double> &output)
	{
		auto dmap = decomposition->dmap->cbegin();
		for (size_t i = 0; i < output.size(); ++i, ++dmap) {
			for (auto di = dmap->begin(); di != dmap->end(); ++di) {
				if (decomposition->ismy(di->domain)) {
					output[i] = domains[di->domain - decomposition->dbegin].vals[di->index];
					break; // we assume synchronization inside the solver
				}
			}
		}
	}

public:
	void storeTo(std::vector<double> &output)
	{
		_store(this->domains, this->decomposition, output);
	}

	void copyReal(const Vector_Distributed<Vector_Dense , std::complex<T> > *a)
	{
		eslog::error("call empty function\n");
	}

	void copyReal(const Vector_Distributed<Vector_Sparse, std::complex<T> > *a)
	{
		eslog::error("call empty function\n");
	}

	void copyReal(const Vector_FETI<Vector_Dense , std::complex<T> > *a)
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(this->domains[d], a->domains[d], 0);
		}
	}

	void copyReal(const Vector_FETI<Vector_Sparse, std::complex<T> > *a)
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(this->domains[d], a->domains[d], 0);
		}
	}

	void copyImag(const Vector_Distributed<Vector_Dense , std::complex<T> > *a)
	{
		eslog::error("call empty function\n");
	}

	void copyImag(const Vector_Distributed<Vector_Sparse, std::complex<T> > *a)
	{
		eslog::error("call empty function\n");
	}

	void copyImag(const Vector_FETI<Vector_Dense , std::complex<T> > *a)
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(this->domains[d], a->domains[d], 1);
		}
	}

	void copyImag(const Vector_FETI<Vector_Sparse, std::complex<T> > *a)
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(this->domains[d], a->domains[d], 1);
		}
	}

	void copyToReal(Vector_Distributed<Vector_Dense , std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToReal(Vector_Distributed<Vector_Sparse, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToReal(Vector_FETI<Vector_Dense , std::complex<T> > *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], 0, this->domains[d]);
		}
	}

	void copyToReal(Vector_FETI<Vector_Sparse, std::complex<T> > *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], 0, this->domains[d]);
		}
	}

	void copyToImag(Vector_Distributed<Vector_Dense , std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToImag(Vector_Distributed<Vector_Sparse, std::complex<T> > *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyToImag(Vector_FETI<Vector_Dense , std::complex<T> > *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], 1, this->domains[d]);
		}
	}

	void copyToImag(Vector_FETI<Vector_Sparse, std::complex<T> > *a) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], 1, this->domains[d]);
		}
	}

	void copySliced(const Vector_Base<T> *in, int offset, int size, int step)
	{
		in->copyToSliced(this, offset, size, step);
	}

	void addSliced(const T &alpha, const Vector_Base<double> *a, int offset, int size, int step)
	{
		a->addToSliced(alpha, this, offset, size, step);
	}

	void copyToSliced(Vector_Distributed<Vector_Dense, T> *a, int offset, int size, int step) const
	{
		eslog::error("call empty function\n");
	}

	void copyToSliced(Vector_Distributed<Vector_Sparse, T> *a, int offset, int size, int step) const
	{
		eslog::error("call empty function\n");
	}

	void copyToSliced(Vector_FETI<Vector_Dense, T> *a, int offset, int size, int step) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], this->domains[d], offset, size, step);
		}
	}

	void copyToSliced(Vector_FETI<Vector_Sparse, T> *a, int offset, int size, int step) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::copy(a->domains[d], this->domains[d], offset, size, step);
		}
	}

	void addToSliced(const T &alpha, Vector_Distributed<Vector_Dense, double> *a, int offset, int size, int step) const
	{
		eslog::error("call empty function\n");
	}

	void addToSliced(const T &alpha, Vector_Distributed<Vector_Sparse, double> *a, int offset, int size, int step) const
	{
		eslog::error("call empty function\n");
	}

	void addToSliced(const T &alpha, Vector_FETI<Vector_Dense, double> *a, int offset, int size, int step) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::add(a->domains[d], alpha, this->domains[d], offset, size, step);
		}
	}

	void addToSliced(const T &alpha, Vector_FETI<Vector_Sparse, double> *a, int offset, int size, int step) const
	{
		#pragma omp parallel for
		for (size_t d = 0; d < this->domains.size(); ++d) {
			math::add(a->domains[d], alpha, this->domains[d], offset, size, step);
		}
	}
};

template <template<typename> typename Vector, typename T>
class Vector_FETI<Vector, std::complex<T> >: public Vector_FETI_Common<Vector, std::complex<T> > {
public:
	void storeTo(std::vector<double> &output)
	{
//		for (size_t i = 0; i < output.size(); ++i) {
//			output[i] = this->cluster.vals[i].real();
//		}
		eslog::error("call empty function\n");
	}

	void copyReal(const Vector_Base<T> *in)
	{
		in->copyToReal(this);
	}

	void copyImag(const Vector_Base<T> *in)
	{
		in->copyToImag(this);
	}

	void copyRealTo(Vector_Base<T> *in) const
	{
		in->copyReal(this);
	}

	void copyImagTo(Vector_Base<T> *in) const
	{
		in->copyImag(this);
	}
};
}

#endif /* SRC_MATH2_GENERALIZATION_VECTOR_FETI_H_ */
