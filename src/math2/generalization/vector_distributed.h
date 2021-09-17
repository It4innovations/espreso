
#ifndef SRC_MATH2_GENERALIZATION_VECTOR_DISTRIBUTED_H_
#define SRC_MATH2_GENERALIZATION_VECTOR_DISTRIBUTED_H_

#include "vector_base.h"
#include "esinfo/eslog.h"
#include "math2/math2.h"
#include "math2/utils/dofs_distribution.h"
#include "wrappers/mpi/communication.h"

#include <vector>

namespace espreso {

template <typename Parent, typename V, template<typename> typename Vector, typename T>
class Vector_Distributed_Common: public Parent {
public:
	Parent* copyPattern()
	{
		V *m = new V();
		m->cluster.pattern(cluster);
		return m;
	}

	void store(const char *file)
	{
		math::store(*static_cast<V*>(this), file);
	}

	void set(const T &value)
	{
		math::set(cluster, value);
	}

	void scale(const T &alpha)
	{
		math::scale(alpha, cluster);
	}

	void copy(const Vector_Base<T> *in)
	{
		if (dynamic_cast<const Vector_Distributed<Vector, T>*>(in)) {
			math::copy(cluster, static_cast<const Vector_Distributed<Vector, T>*>(in)->cluster);
		}
	}

	void add(const T &alpha, const Vector_Base<T> *a)
	{
		a->addTo(alpha, static_cast<V*>(this));
	}

	T norm()
	{
		T dot = math::dot(cluster.size - distribution.halo.size(), cluster.vals + distribution.halo.size(), 1, cluster.vals + distribution.halo.size(), 1);
		Communication::allReduce(&dot, NULL, 1, MPI_DOUBLE, MPI_SUM);
		return std::sqrt(dot);
	}

	T max()
	{
		eslog::error("call empty function: max\n");
		return 0;
	}

	T absmax()
	{
		eslog::error("call empty function: absmax\n");
		return 0;
	}

	T dot(const Vector_Base<T> *other)
	{
		eslog::error("call empty function: dot\n");
		return 0;
	}

	void addTo(const T &alpha, Vector_Distributed<Vector_Dense, T> *a) const
	{
		math::add(a->cluster, alpha, cluster);
	}

	void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, T> *a) const
	{
		math::add(a->cluster, alpha, cluster);
	}

	Vector<T> cluster;
	DOFsDistribution distribution;
};


template <template<typename> typename Vector, typename T>
class Vector_Distributed: public Vector_Distributed_Common<Vector_Base<T>, Vector_Distributed<Vector, T>, Vector, T> {
public:
	void store(std::vector<double> &output)
	{
		for (size_t i = 0; i < output.size(); ++i) {
			output[i] = this->cluster.vals[i];
		}
	}

	void copy(const Vector_Base<T> *in, int offset, int size, int step)
	{
		if (dynamic_cast<const Vector_Distributed<Vector, T>*>(in)) {
			math::copy(this->cluster, static_cast<const Vector_Distributed<Vector, T>*>(in)->cluster, offset, size, step);
		}
	}

	void add(const T &alpha, const Vector_Base<double> *a, int offset, int size, int step)
	{
		a->addTo(alpha, this, offset, size, step);
	}

	void addTo(const T &alpha, Vector_Distributed<Vector_Dense, double> *a, int offset, int size, int step) const
	{
		math::add(a->cluster, alpha, this->cluster, offset, size, step);
	}

	void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, double> *a, int offset, int size, int step) const
	{
		math::add(a->cluster, alpha, this->cluster, offset, size, step);
	}
};

template <template<typename> typename Vector, typename T>
class Vector_Distributed<Vector, std::complex<T> >: public Vector_Distributed_Common<Vector_Base<std::complex<T> >, Vector_Distributed<Vector, std::complex<T> >, Vector, std::complex<T> > {
public:
	void store(std::vector<double> &output)
	{
		for (size_t i = 0; i < output.size(); ++i) {
			output[i] = this->cluster.vals[i].real();
		}
	}

	void copyReal(const Vector_Base<T> *in)
	{
		if (dynamic_cast<const Vector_Distributed<Vector, T >*>(in)) {
			math::copy(this->cluster, 0, static_cast<const Vector_Distributed<Vector, T >*>(in)->cluster);
		}
	}

	void copyImag(const Vector_Base<T> *in)
	{
		if (dynamic_cast<const Vector_Distributed<Vector, T >*>(in)) {
			math::copy(this->cluster, 1, static_cast<const Vector_Distributed<Vector, T >*>(in)->cluster);
		}
	}

	void copyRealTo(Vector_Base<T> *in) const
	{
		if (dynamic_cast<const Vector_Distributed<Vector, T >*>(in)) {
			math::copy(static_cast<Vector_Distributed<Vector, T >*>(in)->cluster, this->cluster, 0);
		}
	}

	void copyImagTo(Vector_Base<T> *in) const
	{
		if (dynamic_cast<const Vector_Distributed<Vector, T >*>(in)) {
			math::copy(static_cast<Vector_Distributed<Vector, T >*>(in)->cluster, this->cluster, 1);
		}
	}
};

}

#endif /* SRC_MATH2_GENERALIZATION_VECTOR_DISTRIBUTED_H_ */
