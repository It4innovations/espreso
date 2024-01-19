
#ifndef SRC_MATH2_GENERALIZATION_VECTOR_FETI_H_
#define SRC_MATH2_GENERALIZATION_VECTOR_FETI_H_

#include "vector_base.h"
#include "matrix_feti.decomposition.h"
#include "math.physics.h"
#include "esinfo/eslog.h"

#include <vector>

namespace espreso {

template <template<typename, typename, template<typename> typename> typename Vector, typename T>
class Vector_FETI: public Vector_Base<T> {

	void _store(double &out, const double &in) { out = in; }
	void _store(double &out, const std::complex<double> &in) { out = in.real(); }

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

	void storeTo(std::vector<double> &output)
	{
		auto dmap = decomposition->dmap->cbegin();
		for (size_t i = 0; i < output.size(); ++i, ++dmap) {
			for (auto di = dmap->begin(); di != dmap->end(); ++di) {
				if (decomposition->ismy(di->domain)) {
					_store(output[i], domains[di->domain - decomposition->dbegin].vals[di->index]);
					break; // we assume synchronization inside the solver
				}
			}
		}
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

	std::vector<Vector<T, int, cpu_allocator> > domains;
	DOFsDecomposition *decomposition;
};

}

#endif /* SRC_MATH2_GENERALIZATION_VECTOR_FETI_H_ */
