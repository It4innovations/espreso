
#include "dual_map.h"
#include "feti/feti.h"

#include <complex>
#include <vector>

namespace espreso {

template <typename T>
void Dual_Map::set(FETI<T> &feti) {
    nhalo = feti.lambdas.nhalo;
    size = feti.lambdas.size;
    neighbors = feti.decomposition->neighbors;

    nsize.resize(neighbors.size());
    for (size_t i = 0, offset = 0; i < feti.lambdas.cmap.size(); ) {
        esint lambdas =  feti.lambdas.cmap[i];
        esint domains =  feti.lambdas.cmap[i + 1];
        nmap.push_back(offset);
        nmap.push_back(offset + lambdas);
        size_t ncounter = nmap.size();
        nmap.push_back(0); // neighbors
        esint last = -1;
        for (esint d = 0; d < domains; ++d) {
            if (!feti.decomposition->ismy(feti.lambdas.cmap[2 + i + d])) {
                esint neigh = feti.decomposition->noffset(feti.lambdas.cmap[2 + i + d]);
                if (last < neigh) {
                    nmap.push_back(neigh);
                    nmap[ncounter]++;
                    nsize[neigh] += lambdas;
                    last = neigh;
                }
            }
        }
        if (nmap[ncounter] == 0) {
            nmap.resize(nmap.size() - 3);
        }
        offset += lambdas;
        i += feti.lambdas.cmap[i + 1] + 2;
    }
}

int Dual_Map::nhalo, Dual_Map::size;
std::vector<int> Dual_Map::nmap, Dual_Map::neighbors, Dual_Map::nsize;

template void Dual_Map::set(FETI<double>&);
template void Dual_Map::set(FETI<std::complex<double> >&);

}

