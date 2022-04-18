
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

	static void set(esint nhalo, esint size)
	{
		Vector_Dual<T>::nhalo = nhalo;
		Vector_Dual<T>::halo.resize(size + align);
		void* _vals = static_cast<void*>(Vector_Dual<T>::halo.vals);
		size_t _size = size + align;
		Vector_Dual<T>::halo.vals = static_cast<T*>(std::align(alignof(T), sizeof(T), _vals, _size));
		Vector_Dual<T>::halo.size = size;
	}

	void resize(esint size)
	{
		Vector_Dense<T>::resize(size + align);
		void* _vals = static_cast<void*>(Vector_Dense<T>::vals);
		size_t _size = size + align;
		Vector_Dense<T>::vals = static_cast<T*>(std::align(alignof(T), sizeof(T), _vals, _size));
		Vector_Dense<T>::size = size;
	}

	void synchronize()
	{
//		Communication::allGatherInplace(this->vals, Vector_Kernel<T>::offset, Vector_Kernel<T>::size);
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
	static esint nhalo;
	static Vector_Dense<T> halo;
};

template <typename T>
esint Vector_Dual<T>::nhalo = 0;
template <typename T>
Vector_Dense<T> Vector_Dual<T>::halo;

}

#endif /* SRC_MATH2_FETI_VECTOR_DUAL_H_ */
