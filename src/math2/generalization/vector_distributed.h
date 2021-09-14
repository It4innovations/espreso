
#ifndef SRC_MATH2_GENERALIZATION_VECTOR_DISTRIBUTED_H_
#define SRC_MATH2_GENERALIZATION_VECTOR_DISTRIBUTED_H_

#include "vector_base.h"
#include "math2/math2.h"
#include "math2/utils/dofs_distribution.h"
#include "wrappers/mpi/communication.h"

#include <vector>

namespace espreso {

template <template<typename> typename Vector, typename T>
class Vector_Distributed: public Vector_Base<T>
{
public:
	~Vector_Distributed()
	{

	}

//	Vector_Distributed* copy()
//	{
//		Vector_Distributed<Vector, T> *m = new Vector_Distributed<Vector, T>();
//		m->cluster.resize(cluster);
//		return m;
//	}

	Vector_Distributed* copyPattern()
	{
		Vector_Distributed<Vector, T> *m = new Vector_Distributed<Vector, T>();
		m->cluster.pattern(cluster);
		return m;
	}

	void store(const char *file)
	{
		math::store(*this, file);
	}

	void store(std::vector<T> &output)
	{
		for (size_t i = 0; i < output.size(); ++i) {
			output[i] = cluster.vals[i];
		}
	}

	void fill(const T &value)
	{
		math::fill(cluster, value);
	}

	void fillData(const Vector_Base<T> *in)
	{
		if (dynamic_cast<const Vector_Distributed<Vector, T>*>(in)) {
			math::copy(cluster, static_cast<const Vector_Distributed<Vector, T>*>(in)->cluster);
		}
	}

	void fillData(const Vector_Base<T> *in, int offset, int size, int step)
	{
		if (dynamic_cast<const Vector_Distributed<Vector, T>*>(in)) {
			math::copy(cluster, static_cast<const Vector_Distributed<Vector, T>*>(in)->cluster, offset, size, step);
		}
	}

	void scale(const T &alpha)
	{
		math::scale(alpha, cluster);
	}

	void add(const T &alpha, const Vector_Base<T> *a)
	{
		if (dynamic_cast<const Vector_Distributed<Vector, T>*>(a)) {
			math::add<T>(cluster, alpha, static_cast<const Vector_Distributed<Vector, T>*>(a)->cluster);
		}
	}

	void add(const T &alpha, const Vector_Base<T> *a, int offset, int size, int step)
	{
		if (dynamic_cast<const Vector_Distributed<Vector, T>*>(a)) {
			math::add<T>(cluster, alpha, static_cast<const Vector_Distributed<Vector, T>*>(a)->cluster, offset, size, step);
		}
	}

	void sum(const T &alpha, const Vector_Base<T> *a, const T &beta, const Vector_Base<T> *b)
	{
		if (dynamic_cast<const Vector_Distributed<Vector, T>*>(a) && dynamic_cast<const Vector_Distributed<Vector, T>*>(b)) {
			math::sum(cluster, alpha, static_cast<const Vector_Distributed<Vector, T>*>(a)->cluster, beta, static_cast<const Vector_Distributed<Vector, T>*>(b)->cluster);
		}
	}

	void sum(const T &alpha, const Vector_Base<T> *a, const T &beta, const Vector_Base<T> *b, int offset, int size, int step)
	{
		if (dynamic_cast<const Vector_Distributed<Vector, T>*>(a) && dynamic_cast<const Vector_Distributed<Vector, T>*>(b)) {
			math::sum(cluster, alpha, static_cast<const Vector_Distributed<Vector, T>*>(a)->cluster, beta, static_cast<const Vector_Distributed<Vector, T>*>(b)->cluster, offset, size, step);
		}
	}

	void addTo(const T &alpha, Vector_Sparse<T> *a) const
	{
		math::add(*a, alpha, cluster);
	}

	T norm()
	{
		T dot = math::dot(cluster.size - distribution.halo.size(), cluster.vals + distribution.halo.size(), 1, cluster.vals + distribution.halo.size(), 1);
		Communication::allReduce(&dot, NULL, 1, MPI_DOUBLE, MPI_SUM);
		return std::sqrt(dot);
	}

	T max()
	{
		return 0;
	}

	T absmax()
	{
		return 0;
	}

	T dot(const Vector_Base<T> *other)
	{
		return 0;
	}

	Vector<T> cluster;
	DOFsDistribution distribution;
};

}

#endif /* SRC_MATH2_GENERALIZATION_VECTOR_DISTRIBUTED_H_ */
