
#include "matrix_dual_sparse.h"
#include "math/wrappers/math.blas.h"
#include "wrappers/mpi/communication.h"

#include <complex>

namespace espreso {

template <typename T>
Matrix_Dual_Sparse<T>::Matrix_Dual_Sparse(const std::vector<int> &nonzeros, const std::vector<int> &distributed)
: Matrix_Dual<T>(nonzeros.size()), nonzeros(nonzeros), distributed(distributed)
{
    std::vector<std::vector<int> > send(Dual_Map::nsize.size(), distributed);
    recvNonzeros.resize(Dual_Map::nsize.size());
    Communication::exchangeUnknownSize(send, recvNonzeros, Dual_Map::neighbors);
}

template <typename T>
void Matrix_Dual_Sparse<T>::synchronize()
{
//    Matrix_Dual<T>::synchronize();
    std::vector<std::vector<T> > send(Dual_Map::nsize.size()), recv(Dual_Map::nsize.size());
    for (size_t i = 0; i < Dual_Map::nsize.size(); ++i) {
        send[i].resize(Dual_Map::nsize[i] * distributed.size());
    }
    std::vector<int> offset(send.size());
    std::fill(offset.begin(), offset.end(), 0);
    for (size_t rn = 0, rd = 0; rn < nonzeros.size(); ++rn) {
        if (rd < distributed.size() && nonzeros[rn] == distributed[rd]) {
            for (size_t i = 0; i < Dual_Map::nmap.size();) {
                for (int n = 0; n < Dual_Map::nmap[i + 2]; ++n) {
                    int ni = Dual_Map::nmap[i + 3 + n];
                    std::copy(vals + rn * ncols + Dual_Map::nmap[i], vals + rn * ncols + Dual_Map::nmap[i + 1], send[ni].data() + offset[ni]);
                    offset[ni] += Dual_Map::nmap[i + 1] - Dual_Map::nmap[i];
                }
                i += Dual_Map::nmap[i + 2] + 3;
            }
            ++rd;
        }
    }
    Communication::exchangeUnknownSize(send, recv, Dual_Map::neighbors);
    std::fill(offset.begin(), offset.end(), 0);
    for (size_t n = 0; n < recv.size(); ++n) {
        for (size_t rn = 0, rd = 0; rn < nonzeros.size(); ++rn) {
            while (rd < recvNonzeros[n].size() && recvNonzeros[n][rd] < nonzeros[rn]) {
                offset[n] += Dual_Map::nsize[n];
                ++rd;
            }
            if (rd < recvNonzeros[n].size() && recvNonzeros[n][rd] == nonzeros[rn]) {
                for (size_t i = 0; i < Dual_Map::nmap.size();) {
                    if (std::any_of(Dual_Map::nmap.data() + i + 3, Dual_Map::nmap.data() + i + 3 + Dual_Map::nmap[i + 2], [&n] (size_t ni) { return ni == n;})) {
                        math::blas::add<T>(Dual_Map::nmap[i + 1] - Dual_Map::nmap[i], vals + rn * ncols + Dual_Map::nmap[i], 1, 1, recv[n].data() + offset[n], 1);
                        offset[n] += Dual_Map::nmap[i + 1] - Dual_Map::nmap[i];
                    }
                    i += Dual_Map::nmap[i + 2] + 3;
                }
                ++rd;
            }
        }
    }
}

template struct Matrix_Dual_Sparse<double>;
template struct Matrix_Dual_Sparse<std::complex<double> >;

}
