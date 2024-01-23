
#ifndef SRC_MATH2_GENERALIZATION_VECTOR_DISTRIBUTED_H_
#define SRC_MATH2_GENERALIZATION_VECTOR_DISTRIBUTED_H_

#include "vector_base.h"
#include "esinfo/eslog.h"
#include "math.physics.h"
#include "matrix_distributed.distribution.h"
#include "matrix_distributed.synchronization.h"
#include "wrappers/mpi/communication.h"

#include <vector>

namespace espreso {

template <template<typename, typename, typename> typename Vector, typename T>
class Vector_Distributed: public Vector_Base<T> {
public:
	void synchronize()
	{
		synchronization->synchronize(*static_cast<Vector_Distributed<Vector, T>*>(this));
	}

	void scatter()
	{
		synchronization->scatterToUpper(*static_cast<Vector_Distributed<Vector, T>*>(this));
	}

	Vector_Base<T>* copyPattern()
	{
		Vector_Distributed<Vector, T> *m = new Vector_Distributed<Vector, T>();
		m->cluster.pattern(cluster);
		m->distribution = this->distribution;
		m->synchronization = this->synchronization;
		return m;
	}

	void store(const char *file)
	{
		math::store(*static_cast<Vector_Distributed<Vector, T>*>(this), file);
	}

	void storeTo(std::vector<double> &output)
	{
		for (size_t i = 0; i < output.size(); ++i) {
			output[i] = this->cluster.vals[i];
		}
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
		this->touched = true;
		in->copyTo(static_cast<Vector_Distributed<Vector, T>*>(this));
	}

	void add(const T &alpha, const Vector_Base<T> *a)
	{
		a->addTo(alpha, static_cast<Vector_Distributed<Vector, T>*>(this));
	}

	T norm()
	{
		T dot = math::blas::dot(cluster.size - distribution->halo.size(), cluster.vals + distribution->halo.size(), 1, cluster.vals + distribution->halo.size(), 1);
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

	void copyTo(Vector_Distributed<Vector_Dense , T> *a) const
	{
		math::copy(a->cluster, cluster);
	}

	void copyTo(Vector_Distributed<Vector_Sparse, T> *a) const
	{
		math::copy(a->cluster, cluster);
	}

	void copyTo(Vector_FETI<Vector_Dense , T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void copyTo(Vector_FETI<Vector_Sparse, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void addTo(const T &alpha, Vector_Distributed<Vector_Dense, T> *a) const
	{
		math::add(a->cluster, alpha, cluster);
	}

	void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, T> *a) const
	{
		math::add(a->cluster, alpha, cluster);
	}

	void addTo(const T &alpha, Vector_FETI<Vector_Dense, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	void addTo(const T &alpha, Vector_FETI<Vector_Sparse, T> *a) const
	{
		eslog::error("call empty function\n");
	}

	Vector<T, esint, cpu_allocator> cluster;
	DOFsDistribution *distribution;
	Data_Synchronization<Vector, T> *synchronization;
};

}

#endif /* SRC_MATH2_GENERALIZATION_VECTOR_DISTRIBUTED_H_ */
