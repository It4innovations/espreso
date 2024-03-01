
#include "vector_dual.h"

#include "math/wrappers/math.blas.h"
#include "wrappers/mpi/communication.h"

#include <complex>
#include <memory>

namespace espreso {

template <typename T>
Vector_Dual<T>::Vector_Dual()
{
    Vector_Dense<T>::resize(Dual_Map::size);
}

template <typename T>
void Vector_Dual<T>::synchronize()
{
    std::vector<esint> offset(send.size());
    for (size_t i = 0; i < Dual_Map::nmap.size();) {
        for (esint n = 0; n < Dual_Map::nmap[i + 2]; ++n) {
            esint ni = Dual_Map::nmap[i + 3 + n];
            std::copy(vals + Dual_Map::nmap[i], vals + Dual_Map::nmap[i + 1], send[ni].data() + offset[ni]);
            offset[ni] += Dual_Map::nmap[i + 1] - Dual_Map::nmap[i];
        }
        i += Dual_Map::nmap[i + 2] + 3;
    }
    Communication::exchangeKnownSize(send, recv, Dual_Map::neighbors);
    std::fill(offset.begin(), offset.end(), 0);
    for (size_t i = 0; i < Dual_Map::nmap.size();) {
        for (esint n = 0; n < Dual_Map::nmap[i + 2]; ++n) {
            esint ni = Dual_Map::nmap[i + 3 + n];
            math::blas::add<T>(Dual_Map::nmap[i + 1] - Dual_Map::nmap[i], vals + Dual_Map::nmap[i], 1, 1, recv[ni].data() + offset[ni], 1);
            offset[ni] += Dual_Map::nmap[i + 1] - Dual_Map::nmap[i];
        }
        i += Dual_Map::nmap[i + 2] + 3;
    }
}

template <typename T>
void Vector_Dual<T>::copyToWithoutHalo(Vector_Dense<T> &to) const
{
    const auto &li = Dual_Map::local_intervals;
    for (size_t i = 0, prev = 0; i < li.size(); ++i) {
        std::fill(to.vals + prev, to.vals + li[i].start, 0);
        math::blas::copy(li[i].size, to.vals + li[i].start, 1, vals + li[i].start, 1);
        prev = li[i].end;
    }
}

template <typename T>
T Vector_Dual<T>::dot(const Vector_Dense<T> &other) const
{
    T sum = T{0};
    const auto &li = Dual_Map::local_intervals;
    for (size_t i = 0; i < li.size(); ++i) {
        sum += math::blas::dot(li[i].size, vals + li[i].start, 1, other.vals + li[i].start, 1);
    }
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
