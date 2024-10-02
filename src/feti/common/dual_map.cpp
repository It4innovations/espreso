
#include "dual_map.h"
#include "feti/feti.h"
#include "wrappers/mpi/communication.h"

#include <complex>
#include <vector>

namespace espreso {

template <typename T>
void Dual_Map::set(FETI<T> &feti) {
    size = feti.lambdas.size;
    Communication::allReduce(&size, &total, 1, MPITools::getType(size).mpitype, MPI_SUM);
    local_intervals.resize(feti.lambdas.intervals.size());
    for (size_t i = 0, offset = 0; i < local_intervals.size(); ++i) {
        local_intervals[i].start = feti.lambdas.intervals[i].halo + offset;
        local_intervals[i].end = local_intervals[i].start + feti.lambdas.intervals[i].size;
        local_intervals[i].size = local_intervals[i].end - local_intervals[i].start;
        local_intervals[i].offset = local_intervals[i].size;
        local_intervals[i].total = Communication::exscan(local_intervals[i].offset);
        offset = local_intervals[i].end;
    }

    neighbors = feti.decomposition->neighbors;

    nsize.clear();
    nsize.resize(neighbors.size());
    nmap.clear();
    for (size_t i = 0, offset = 0; i < feti.lambdas.cmap.size(); ) {
        int lambdas =  feti.lambdas.cmap[i];
        int domains =  feti.lambdas.cmap[i + 1];
        nmap.push_back(offset);
        nmap.push_back(offset + lambdas);
        size_t ncounter = nmap.size();
        nmap.push_back(0); // neighbors
        int last = -1;
        for (int d = 0; d < domains; ++d) {
            if (!feti.decomposition->ismy(feti.lambdas.cmap[2 + i + d])) {
                int neigh = feti.decomposition->noffset(feti.lambdas.cmap[2 + i + d]);
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

int Dual_Map::size;
int Dual_Map::total;
std::vector<Dual_Map::interval> Dual_Map::local_intervals;
std::vector<int> Dual_Map::nmap, Dual_Map::neighbors, Dual_Map::nsize;

template void Dual_Map::set(FETI<double>&);
template void Dual_Map::set(FETI<std::complex<double> >&);

}

