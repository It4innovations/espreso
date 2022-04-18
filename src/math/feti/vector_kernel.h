
#ifndef SRC_MATH2_FETI_VECTOR_KERNEL_H_
#define SRC_MATH2_FETI_VECTOR_KERNEL_H_

#include "esinfo/envinfo.h"
#include "math/primitives/vector_dense.h"
#include "wrappers/mpi/communication.h"

#include <vector>
#include <memory>

namespace espreso {

template <typename T>
struct Vector_Kernel: public Vector_Dense<T> {

	enum: size_t {
		align = 64U
	};

	static void set(esint offset, esint size, esint totalSize)
	{
		Vector_Kernel<T>::offset = offset;
		Vector_Kernel<T>::size = size;
		Vector_Kernel<T>::totalSize = totalSize;
		Vector_Kernel<T>::distribution.resize(info::env::threads + 1);
		size_t chunk = align / sizeof(T);
		size_t tsize = size / chunk;
		for (size_t t = 1; t < Vector_Kernel<T>::distribution.size(); ++t) {
			Vector_Kernel<T>::distribution[t] = Vector_Kernel<T>::distribution[t - 1] + tsize * chunk;
			if (size % chunk < t - 1) {
				Vector_Kernel<T>::distribution[t] += chunk;
			}
		}
		Vector_Kernel<T>::distribution.back() = size;
	}

	void resize()
	{
		Vector_Dense<T>::resize(Vector_Kernel<T>::totalSize + 2 * align);
		void* _vals = static_cast<void*>(Vector_Dense<T>::vals);
		size_t _size = Vector_Kernel<T>::totalSize + align;
		Vector_Dense<T>::vals = static_cast<T*>(std::align(align, sizeof(T), _vals, _size));
		Vector_Dense<T>::size = Vector_Kernel<T>::totalSize;
	}

	void synchronize()
	{
		Communication::allGatherInplace(this->vals, Vector_Kernel<T>::offset, Vector_Kernel<T>::size);
	}

	static esint offset, size, totalSize;
	static std::vector<size_t> distribution;
};

template <typename T>
esint Vector_Kernel<T>::offset = 0;
template <typename T>
esint Vector_Kernel<T>::size = 0;
template <typename T>
esint Vector_Kernel<T>::totalSize = 0;
template <typename T>
std::vector<size_t> Vector_Kernel<T>::distribution = { 0, 0 };

}

#endif /* SRC_MATH2_FETI_VECTOR_KERNEL_H_ */
