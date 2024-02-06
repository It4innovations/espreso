
#include "vector_dual.h"

#include "esinfo/envinfo.h"
#include "analysis/builder/feti.decomposition.h"
#include "math/primitives/vector_dense.h"
#include "math/wrappers/math.blas.h"
#include "wrappers/mpi/communication.h"

#include <complex>
#include <memory>

namespace espreso {

template <typename T>
void Vector_Dual<T>::set(esint nhalo, const std::vector<esint> &cmap, const FETIDecomposition &decomposition)
{
	Vector_Dual<T>::nhalo = nhalo;
	Vector_Dual<T>::neighbors.assign(decomposition.neighbors.begin(), decomposition.neighbors.end());
	Vector_Dual<T>::sBuffer.resize(Vector_Dual<T>::neighbors.size());
	Vector_Dual<T>::rBuffer.resize(Vector_Dual<T>::neighbors.size());

	std::vector<esint> size(Vector_Dual<T>::neighbors.size());
	for (size_t i = 0; i < cmap.size(); ) {
		esint lambdas = cmap[i];
		esint domains = cmap[i + 1];
		Vector_Dual<T>::nmap.push_back(Vector_Dual<T>::localSize);
		Vector_Dual<T>::nmap.push_back(Vector_Dual<T>::localSize + lambdas);
		size_t ncounter = Vector_Dual<T>::nmap.size();
		Vector_Dual<T>::nmap.push_back(0); // neighbors
		esint last = -1;
		for (esint d = 0; d < domains; ++d) {
			if (!decomposition.ismy(cmap[2 + i + d])) {
				esint neigh = decomposition.noffset(cmap[2 + i + d]);
				if (last < neigh) {
					Vector_Dual<T>::nmap.push_back(neigh);
					Vector_Dual<T>::nmap[ncounter]++;
					size[neigh] += lambdas;
					last = neigh;
				}
			}
		}
		if (Vector_Dual<T>::nmap[ncounter] == 0) {
			Vector_Dual<T>::nmap.resize(Vector_Dual<T>::nmap.size() - 3);
		}
		Vector_Dual<T>::localSize += lambdas;
		i += cmap[i + 1] + 2;
	}

	for (size_t i = 0; i < Vector_Dual<T>::neighbors.size(); ++i) {
		Vector_Dual<T>::sBuffer[i].resize(size[i]);
		Vector_Dual<T>::rBuffer[i].resize(size[i]);
	}
}

template <typename T>
void Vector_Dual<T>::resize()
{
	Vector_Dense<T>::resize(Vector_Dual<T>::localSize);
}

template <typename T>
void Vector_Dual<T>::synchronize()
{
	std::vector<esint> offset(sBuffer.size());
	for (size_t i = 0; i < nmap.size();) {
		for (esint n = 0; n < nmap[i + 2]; ++n) {
			esint ni = nmap[i + 3 + n];
			std::copy(this->vals + nmap[i], this->vals + nmap[i + 1], sBuffer[ni].data() + offset[ni]);
			offset[ni] += nmap[i + 1] - nmap[i];
		}
		i += nmap[i + 2] + 3;
	}
	Communication::exchangeKnownSize(sBuffer, rBuffer, neighbors);
	std::fill(offset.begin(), offset.end(), 0);
	for (size_t i = 0; i < nmap.size();) {
		for (esint n = 0; n < nmap[i + 2]; ++n) {
			esint ni = nmap[i + 3 + n];
			math::blas::add<T>(nmap[i + 1] - nmap[i], this->vals + nmap[i], 1, 1, rBuffer[ni].data() + offset[ni], 1);
			offset[ni] += nmap[i + 1] - nmap[i];
		}
		i += nmap[i + 2] + 3;
	}
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
		to.vals[i] = 0;
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
float Vector_Dual<float>::dot(const Vector_Dense<float> &other) const
{
	float sum = 0;
	#pragma omp parallel for reduction(+:sum)
	for (esint i = nhalo; i < this->size; ++i) {
		sum += other.vals[i] * this->vals[i];
	}
	Communication::allReduce(&sum, nullptr, 1, MPITools::getType(sum).mpitype, MPI_SUM);
	return sum;
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
std::complex<float> Vector_Dual<std::complex<float> >::dot(const Vector_Dense<std::complex<float> > &other) const
{
	std::complex<float> sum = 0;
	return sum;
}

template <>
std::complex<double> Vector_Dual<std::complex<double> >::dot(const Vector_Dense<std::complex<double> > &other) const
{
	std::complex<double> sum = 0;
	return sum;
}

template <typename T>
T Vector_Dual<T>::dot() const
{
	return dot(*this);
}

template struct Vector_Dual<float>;
template struct Vector_Dual<double>;
template struct Vector_Dual<std::complex<float>>;
template struct Vector_Dual<std::complex<double>>;

}
