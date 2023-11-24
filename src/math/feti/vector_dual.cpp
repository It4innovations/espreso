
#include "vector_dual.h"

#include "esinfo/envinfo.h"
#include "math/primitives/vector_dense.h"
#include "math/feti/lmap.h"
#include "math/wrappers/math.blas.h"
#include "wrappers/mpi/communication.h"

#include <complex>
#include <memory>

namespace espreso {

template <typename T>
void Vector_Dual<T>::set(esint nhalo, esint size, const std::vector<LMAP> &lmap, const std::vector<int> &neighbors)
{
	Vector_Dual<T>::nhalo = nhalo;
	Vector_Dual<T>::localSize = size;
	Vector_Dual<T>::halo.resize(size + alignof(T));
	Vector_Dual<T>::neighbors = neighbors;
	void* _vals = static_cast<void*>(Vector_Dual<T>::halo.vals);
	size_t _size = size + alignof(T);
	Vector_Dual<T>::halo.vals = static_cast<T*>(std::align(alignof(T), sizeof(T), _vals, _size));
	Vector_Dual<T>::halo.size = size;

	Vector_Dual<T>::nshared = size;
	Vector_Dual<T>::distribution.resize(neighbors.size() + 1, size);
	esint n = 0;
	for (auto map = lmap.begin(); map != lmap.end(); ++map) {
		if (map->neigh == LMAP::LOCAL || map->neigh == LMAP::DIRICHLET) {
			Vector_Dual<T>::nshared = map->offset;
			break;
		}
		Vector_Dual<T>::distribution[map->neigh] = std::min(Vector_Dual<T>::distribution[map->neigh], map->offset);
		n = map->neigh;
	}
	for (size_t nn = n; nn < neighbors.size(); ++nn) {
		Vector_Dual<T>::distribution[nn + 1] = Vector_Dual<T>::nshared;
	}
}

template <typename T>
void Vector_Dual<T>::resize()
{
	// aligned header + space at the end
	Vector_Dense<T>::resize(Vector_Dual<T>::localSize + 2 * alignof(T));
	void* _vals = static_cast<void*>(Vector_Dense<T>::vals);
	size_t _size = Vector_Dense<T>::size;
	Vector_Dense<T>::vals = static_cast<T*>(std::align(alignof(T), sizeof(T), _vals, _size));
	Vector_Dense<T>::size = Vector_Dual<T>::localSize;
}

template <typename T>
void Vector_Dual<T>::synchronize()
{
	Communication::exchangeKnownSizeInPlace(this->vals, Vector_Dual<T>::halo.vals, Vector_Dual<T>::distribution, Vector_Dual<T>::neighbors);
	math::add<T>(Vector_Dual<T>::nshared, this->vals, 1, 1, Vector_Dual<T>::halo.vals, 1);
}

template <typename T>
void Vector_Dual<T>::copyTo(Vector_Dense<T> &to) const
{
	#pragma omp parallel for
	for (esint i = 0; i < this->size; ++i) {
		to.vals[i] = this->vals[i];
	}

	// check performance
//		#pragma omp parallel for
//		for (esint i = 0; i < this->size; ++i) {
//			math::copy(to, *this);
//		}
}

template <typename T>
void Vector_Dual<T>::copyToWithoutHalo(Vector_Dense<T> &to) const
{
	#pragma omp parallel for
	for (esint i = 0; i < nhalo; ++i) {
		to.vals[i] = T{0};
	}
	#pragma omp parallel for
	for (esint i = nhalo; i < this->size; ++i) {
		to.vals[i] = this->vals[i];
	}
}

template <typename T>
void Vector_Dual<T>::scale(const T &alpha)
{
	#pragma omp parallel for
	for (esint i = 0; i < this->size; ++i) {
		this->vals[i] *= alpha;
	}
}

template <typename T>
void Vector_Dual<T>::add(const T &alpha, const Vector_Dense<T> &other)
{
	#pragma omp parallel for
	for (esint i = 0; i < this->size; ++i) {
		this->vals[i] += alpha * other.vals[i];
	}
}

template <>
double Vector_Dual<double>::dot(const Vector_Dense<double> &other) const
{
	double sum = 0;
	#pragma omp parallel for reduction(+:sum)
	for (esint i = nhalo; i < this->size; ++i) {
		sum += other.vals[i] * this->vals[i];
	}
	Communication::allReduce(&sum, nullptr, 1, MPITools::getType(sum).mpitype, MPI_SUM);
	return sum;
}

template <>
double Vector_Dual<std::complex<double> >::dot(const Vector_Dense<std::complex<double> > &other) const
{
	double sum = 0;
//	#pragma omp parallel for reduction(+:sum)
//	for (esint i = nhalo; i < this->size; ++i) {
//		sum += other.vals[i] * this->vals[i];
//	}
//	Communication::allReduce(&sum, nullptr, 1, MPITools::getType(sum).mpitype, MPI_SUM);
	return sum;
}

template <typename T>
double Vector_Dual<T>::dot() const
{
	return dot(*this);
}

template struct Vector_Dual<double>;
template struct Vector_Dual<std::complex<double> >;

}
