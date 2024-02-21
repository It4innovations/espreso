
#include "vector_kernel.h"
#include "math/math.h"

#include <complex>

namespace espreso {

template <typename T>
void Vector_Kernel<T>::set(esint offset, esint size, esint totalSize)
{
    Vector_Kernel<T>::offset = offset;
    Vector_Kernel<T>::localSize = size;
    Vector_Kernel<T>::totalSize = totalSize;
    Vector_Kernel<T>::distribution.resize(info::env::threads + 1);
    size_t chunk = align / sizeof(T);
    size_t tsize = std::max(size / chunk, 1UL);
    for (size_t t = 1; t < Vector_Kernel<T>::distribution.size(); ++t) {
        Vector_Kernel<T>::distribution[t] = Vector_Kernel<T>::distribution[t - 1] + tsize * chunk;
        if (size % chunk < t - 1) {
            Vector_Kernel<T>::distribution[t] += chunk;
        }
        Vector_Kernel<T>::distribution[t] = std::min(Vector_Kernel<T>::distribution[t], (size_t)size);
    }
    Vector_Kernel<T>::distribution.back() = size;
}

template <typename T>
Vector_Kernel<T>::Vector_Kernel()
{
    Vector_Dense<T>::resize(Vector_Kernel<T>::totalSize);
}

template <typename T>
void Vector_Kernel<T>::synchronize()
{
    Communication::allGatherInplace(this->vals, Vector_Kernel<T>::offset, Vector_Kernel<T>::localSize);
}

template <typename T>
T Vector_Kernel<T>::dot(const Vector_Dense<T> &other) const
{
    T sum = math::blas::dot(Vector_Kernel<T>::localSize, vals + Vector_Kernel<T>::offset, 1, other.vals + Vector_Kernel<T>::offset, 1);
    Communication::allReduce(&sum, nullptr, 1, MPITools::getType(sum).mpitype, MPI_SUM);
    return sum;
}

template <typename T>
T Vector_Kernel<T>::dot() const
{
    return dot(*this);
}

template struct Vector_Kernel<double>;
template struct Vector_Kernel<std::complex<double> >;

}
