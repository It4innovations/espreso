
#include "vector_dual.h"

#include "analysis/builder/feti.decomposition.h"
#include "math/primitives/vector_dense.h"
#include "math/wrappers/math.blas.h"
#include "wrappers/mpi/communication.h"

#include <complex>
#include <memory>

namespace espreso {

template <typename T>
Vector_Dual<T>::Vector_Dual()
{
	this->nhalo = Dual_Buffer<T>::nhalo;
	Vector_Dense<T>::resize(Dual_Buffer<T>::size);
}

template <typename T>
void Vector_Dual<T>::synchronize()
{
	std::vector<esint> offset(Dual_Buffer<T>::send.size());
	for (size_t i = 0; i < Dual_Buffer<T>::nmap.size();) {
		for (esint n = 0; n < Dual_Buffer<T>::nmap[i + 2]; ++n) {
			esint ni = Dual_Buffer<T>::nmap[i + 3 + n];
			std::copy(this->vals + Dual_Buffer<T>::nmap[i], this->vals + Dual_Buffer<T>::nmap[i + 1], Dual_Buffer<T>::send[ni].data() + offset[ni]);
			offset[ni] += Dual_Buffer<T>::nmap[i + 1] - Dual_Buffer<T>::nmap[i];
		}
		i += Dual_Buffer<T>::nmap[i + 2] + 3;
	}
	Communication::exchangeKnownSize(Dual_Buffer<T>::send, Dual_Buffer<T>::recv, Dual_Buffer<T>::neighbors);
	std::fill(offset.begin(), offset.end(), 0);
	for (size_t i = 0; i < Dual_Buffer<T>::nmap.size();) {
		for (esint n = 0; n < Dual_Buffer<T>::nmap[i + 2]; ++n) {
			esint ni = Dual_Buffer<T>::nmap[i + 3 + n];
			math::blas::add<T>(Dual_Buffer<T>::nmap[i + 1] - Dual_Buffer<T>::nmap[i], this->vals + Dual_Buffer<T>::nmap[i], 1, 1, Dual_Buffer<T>::recv[ni].data() + offset[ni], 1);
			offset[ni] += Dual_Buffer<T>::nmap[i + 1] - Dual_Buffer<T>::nmap[i];
		}
		i += Dual_Buffer<T>::nmap[i + 2] + 3;
	}
}

template <typename T>
void Vector_Dual<T>::copyToWithoutHalo(Vector_Dense<T> &to) const
{
	std::fill(to.vals, to.vals + nhalo, 0);
	math::blas::copy(size - nhalo, to.vals + nhalo, 1, vals + nhalo, 1);
}

template <typename T>
T Vector_Dual<T>::dot(const Vector_Dense<T> &other) const
{
	T sum = math::blas::dot(size - nhalo, vals + nhalo, 1, other.vals + nhalo, 1);
	Communication::allReduce(&sum, nullptr, 1, MPITools::getType(sum).mpitype, MPI_SUM);
	return sum;
}

template <typename T>
T Vector_Dual<T>::dot() const
{
	return dot(*this);
}

template struct Vector_Dual<int>;
template struct Vector_Dual<double>;
template struct Vector_Dual<std::complex<double> >;

}
