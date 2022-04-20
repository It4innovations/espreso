
#ifndef SRC_MATH2_FETI_VECTOR_DUAL_H_
#define SRC_MATH2_FETI_VECTOR_DUAL_H_

#include "esinfo/envinfo.h"
#include "math/primitives/vector_dense.h"
#include "wrappers/mpi/communication.h"

#include <vector>
#include <memory>

namespace espreso {

template <typename T>
struct Vector_Dual: public Vector_Dense<T> {

	enum: size_t {
		align = 8U
	};

	static void set(esint nhalo, esint size, const std::vector<LMAP> &lmap, const std::vector<int> &neighbors)
	{
		Vector_Dual<T>::nhalo = nhalo;
		Vector_Dual<T>::localSize = size;
		Vector_Dual<T>::halo.resize(size + align);
		Vector_Dual<T>::neighbors = neighbors;
		void* _vals = static_cast<void*>(Vector_Dual<T>::halo.vals);
		size_t _size = size + align;
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

	void resize()
	{
		Vector_Dense<T>::resize(Vector_Dual<T>::localSize + align);
		void* _vals = static_cast<void*>(Vector_Dense<T>::vals);
		size_t _size = Vector_Dual<T>::localSize + align;
		Vector_Dense<T>::vals = static_cast<T*>(std::align(alignof(T), sizeof(T), _vals, _size));
		Vector_Dense<T>::size = Vector_Dual<T>::localSize;
	}

	void synchronize()
	{
		Communication::exchangeKnownSizeInPlace(this->vals, Vector_Dual<T>::halo.vals, Vector_Dual<T>::distribution, Vector_Dual<T>::neighbors);
		math::add<T>(Vector_Dual<T>::nshared, this->vals, 1, 1, Vector_Dual<T>::halo.vals, 1);
	}

	void copyTo(Vector_Dense<T> &to) const
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

	void copyToWithoutHalo(Vector_Dense<T> &to) const
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

	void scale(const T &alpha)
	{
		#pragma omp parallel for
		for (esint i = 0; i < this->size; ++i) {
			this->vals[i] *= alpha;
		}
	}

	void add(const T &alpha, const Vector_Dense<T> &other)
	{
		#pragma omp parallel for
		for (esint i = 0; i < this->size; ++i) {
			this->vals[i] += alpha * other.vals[i];
		}
	}

	double dot(const Vector_Dense<T> &other) const
	{
		double sum = 0;
		#pragma omp parallel for reduction(+:sum)
		for (esint i = nhalo; i < this->size; ++i) {
			sum += other.vals[i] * this->vals[i];
		}
		Communication::allReduce(&sum, nullptr, 1, MPITools::getType(sum).mpitype, MPI_SUM);
		return sum;
	}

	double dot() const
	{
		return dot(*this);
	}

protected:
	static esint nhalo, nshared, localSize;
	static Vector_Dense<T> halo;
	static std::vector<esint> distribution;
	static std::vector<int> neighbors;
};

template <typename T>
esint Vector_Dual<T>::nhalo = 0;
template <typename T>
esint Vector_Dual<T>::nshared = 0;
template <typename T>
esint Vector_Dual<T>::localSize = 0;
template <typename T>
Vector_Dense<T> Vector_Dual<T>::halo;
template <typename T>
std::vector<esint> Vector_Dual<T>::distribution;
template <typename T>
std::vector<esint> Vector_Dual<T>::neighbors;

}

#endif /* SRC_MATH2_FETI_VECTOR_DUAL_H_ */
