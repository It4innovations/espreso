
#ifndef SRC_MATH2_GENERALIZATION_VECTOR_DISTRIBUTED_H_
#define SRC_MATH2_GENERALIZATION_VECTOR_DISTRIBUTED_H_

#include "vector_base.h"
#include "esinfo/eslog.h"
#include "math2/math2.h"
#include "math2/utils/distributed/distribution.h"
#include "math2/utils/distributed/synchronization.h"
#include "wrappers/mpi/communication.h"

#include <vector>

namespace espreso {

template <template<typename> typename Vector, typename T>
class Vector_Distributed_Common: public Vector_Base<T> {
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
		in->copyTo(static_cast<Vector_Distributed<Vector, T>*>(this));
	}

	void add(const T &alpha, const Vector_Base<T> *a)
	{
		a->addTo(alpha, static_cast<Vector_Distributed<Vector, T>*>(this));
	}

	T norm()
	{
		T dot = math::dot(cluster.size - distribution->halo.size(), cluster.vals + distribution->halo.size(), 1, cluster.vals + distribution->halo.size(), 1);
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

	Vector<T> cluster;
	DOFsDistribution *distribution;
	Data_Synchronization<Vector, T> *synchronization;
};


template <template<typename> typename Vector, typename T>
class Vector_Distributed: public Vector_Distributed_Common<Vector, T> {
public:
	void store(std::vector<double> &output)
	{
		for (size_t i = 0; i < output.size(); ++i) {
			output[i] = this->cluster.vals[i];
		}
	}

	void copyReal(const Vector_Distributed<Vector_Dense , std::complex<T> > *a)
	{
		math::copy(this->cluster, a->cluster, 0);
	}

	void copyReal(const Vector_Distributed<Vector_Sparse, std::complex<T> > *a)
	{
		math::copy(this->cluster, a->cluster, 0);
	}

	void copyReal(const Vector_FETI<Vector_Dense , std::complex<T> > *a)
	{

	}

	void copyReal(const Vector_FETI<Vector_Sparse, std::complex<T> > *a)
	{

	}

	void copyImag(const Vector_Distributed<Vector_Dense , std::complex<T> > *a)
	{
		math::copy(this->cluster, a->cluster, 1);
	}

	void copyImag(const Vector_Distributed<Vector_Sparse, std::complex<T> > *a)
	{
		math::copy(this->cluster, a->cluster, 1);
	}

	void copyImag(const Vector_FETI<Vector_Dense , std::complex<T> > *a)
	{

	}

	void copyImag(const Vector_FETI<Vector_Sparse, std::complex<T> > *a)
	{

	}

	void copyToReal(Vector_Distributed<Vector_Dense , std::complex<T> > *a) const
	{
		math::copy(a->cluster, 0, this->cluster);
	}

	void copyToReal(Vector_Distributed<Vector_Sparse, std::complex<T> > *a) const
	{
		math::copy(a->cluster, 0, this->cluster);
	}

	void copyToReal(Vector_FETI<Vector_Dense , std::complex<T> > *a) const
	{

	}

	void copyToReal(Vector_FETI<Vector_Sparse, std::complex<T> > *a) const
	{

	}

	void copyToImag(Vector_Distributed<Vector_Dense , std::complex<T> > *a) const
	{
		math::copy(a->cluster, 1, this->cluster);
	}

	void copyToImag(Vector_Distributed<Vector_Sparse, std::complex<T> > *a) const
	{
		math::copy(a->cluster, 1, this->cluster);
	}

	void copyToImag(Vector_FETI<Vector_Dense , std::complex<T> > *a) const
	{

	}

	void copyToImag(Vector_FETI<Vector_Sparse, std::complex<T> > *a) const
	{

	}

	void copy(const Vector_Base<T> *in, int offset, int size, int step)
	{
		in->copyTo(this, offset, size, step);
	}

	void add(const T &alpha, const Vector_Base<double> *a, int offset, int size, int step)
	{
		a->addTo(alpha, this, offset, size, step);
	}

	void copyTo(Vector_Distributed<Vector_Dense, T> *a, int offset, int size, int step) const
	{
		math::copy(a->cluster, this->cluster, offset, size, step);
	}

	void copyTo(Vector_Distributed<Vector_Sparse, T> *a, int offset, int size, int step) const
	{
		math::copy(a->cluster, this->cluster, offset, size, step);
	}

	void copyTo(Vector_FETI<Vector_Dense, T> *a, int offset, int size, int step) const
	{

	}

	void copyTo(Vector_FETI<Vector_Sparse, T> *a, int offset, int size, int step) const
	{

	}

	void addTo(const T &alpha, Vector_Distributed<Vector_Dense, double> *a, int offset, int size, int step) const
	{
		math::add(a->cluster, alpha, this->cluster, offset, size, step);
	}

	void addTo(const T &alpha, Vector_Distributed<Vector_Sparse, double> *a, int offset, int size, int step) const
	{
		math::add(a->cluster, alpha, this->cluster, offset, size, step);
	}

	void addTo(const T &alpha, Vector_FETI<Vector_Dense, double> *a, int offset, int size, int step) const
	{

	}

	void addTo(const T &alpha, Vector_FETI<Vector_Sparse, double> *a, int offset, int size, int step) const
	{

	}
};

template <template<typename> typename Vector, typename T>
class Vector_Distributed<Vector, std::complex<T> >: public Vector_Distributed_Common<Vector, std::complex<T> > {
public:
	void store(std::vector<double> &output)
	{
		for (size_t i = 0; i < output.size(); ++i) {
			output[i] = this->cluster.vals[i].real();
		}
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

#endif /* SRC_MATH2_GENERALIZATION_VECTOR_DISTRIBUTED_H_ */
