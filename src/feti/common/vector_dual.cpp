
#include "vector_dual.h"

#include "analysis/builder/feti.decomposition.h"
#include "math/primitives/vector_dense.h"
#include "math/wrappers/math.blas.h"
#include "wrappers/mpi/communication.h"

#include <complex>
#include <memory>

namespace espreso {

template <typename T>
void Vector_Dual<T>::resize(int nhalo, int size)
{
	this->nhalo = nhalo;
	Vector_Dense<T>::resize(size);
}

template <typename T>
void Vector_Dual<T>::synchronize(Dual_Buffer<T> &buffer)
{
	std::vector<esint> offset(buffer.send.size());
	for (size_t i = 0; i < buffer.nmap.size();) {
		for (esint n = 0; n < buffer.nmap[i + 2]; ++n) {
			esint ni = buffer.nmap[i + 3 + n];
			std::copy(this->vals + buffer.nmap[i], this->vals + buffer.nmap[i + 1], buffer.send[ni].data() + offset[ni]);
			offset[ni] += buffer.nmap[i + 1] - buffer.nmap[i];
		}
		i += buffer.nmap[i + 2] + 3;
	}
	Communication::exchangeKnownSize(buffer.send, buffer.recv, buffer.neighbors);
	std::fill(offset.begin(), offset.end(), 0);
	for (size_t i = 0; i < buffer.nmap.size();) {
		for (esint n = 0; n < buffer.nmap[i + 2]; ++n) {
			esint ni = buffer.nmap[i + 3 + n];
			math::blas::add<T>(buffer.nmap[i + 1] - buffer.nmap[i], this->vals + buffer.nmap[i], 1, 1, buffer.recv[ni].data() + offset[ni], 1);
			offset[ni] += buffer.nmap[i + 1] - buffer.nmap[i];
		}
		i += buffer.nmap[i + 2] + 3;
	}
}

template <typename T>
void Vector_Dual<T>::copyToWithoutHalo(Vector_Dense<T> &to)
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
template struct Vector_Dual<std::complex<double>>;

}
