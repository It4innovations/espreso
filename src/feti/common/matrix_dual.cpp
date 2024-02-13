
#include "matrix_dual.h"

#include "math/wrappers/math.blas.h"
#include "wrappers/mpi/communication.h"

#include <complex>

namespace espreso {

template <typename T>
void Matrix_Dual<T>::resize(int nrhs)
{
    this->nhalo = Dual_Map::nhalo;
    Matrix_Dense<T>::resize(nrhs, Dual_Map::size);
}

template <typename T>
void Matrix_Dual<T>::synchronize()
{
    if (send.size() == 0 || send[0].size() < (size_t)Dual_Map::nsize[0] * nrows) {
        send.resize(Dual_Map::nsize.size());
        recv.resize(Dual_Map::nsize.size());
        for (size_t i = 0; i < Dual_Map::nsize.size(); ++i) {
            send[i].resize(Dual_Map::nsize[i] * nrows);
            recv[i].resize(Dual_Map::nsize[i] * nrows);
        }
    }
    std::vector<esint> offset(send.size());
    for (int r = 0; r < nrows; ++r) { // is there better than per row solution?
        std::fill(offset.begin(), offset.end(), 0);
        for (size_t i = 0; i < Dual_Map::nmap.size();) {
            for (esint n = 0; n < Dual_Map::nmap[i + 2]; ++n) {
                esint ni = Dual_Map::nmap[i + 3 + n];
                std::copy(vals + r * ncols + Dual_Map::nmap[i], vals + r * ncols + Dual_Map::nmap[i + 1], send[ni].data() + offset[ni]);
                offset[ni] += Dual_Map::nmap[i + 1] - Dual_Map::nmap[i];
            }
            i += Dual_Map::nmap[i + 2] + 3;
        }
    }
    Communication::exchangeKnownSize(send, recv, Dual_Map::neighbors);
    for (int r = 0; r < nrows; ++r) {
        std::fill(offset.begin(), offset.end(), 0);
        for (size_t i = 0; i < Dual_Map::nmap.size();) {
            for (esint n = 0; n < Dual_Map::nmap[i + 2]; ++n) {
                esint ni = Dual_Map::nmap[i + 3 + n];
                math::blas::add<T>(Dual_Map::nmap[i + 1] - Dual_Map::nmap[i], vals + r * ncols + Dual_Map::nmap[i], 1, 1, recv[ni].data() + offset[ni], 1);
                offset[ni] += Dual_Map::nmap[i + 1] - Dual_Map::nmap[i];
            }
            i += Dual_Map::nmap[i + 2] + 3;
        }
    }
}

template struct Matrix_Dual<double>;
template struct Matrix_Dual<std::complex<double> >;

}
