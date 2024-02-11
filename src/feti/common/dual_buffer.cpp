
#include "dual_buffer.h"
#include "feti/feti.h"

namespace espreso {

template <typename T>
void Dual_Buffer<T>::set(FETI<T> &feti)
{
    nhalo = feti.lambdas.nhalo;
    size = feti.lambdas.size;
    neighbors = feti.decomposition->neighbors;
	send.resize(neighbors.size());
	recv.resize(neighbors.size());

	std::vector<esint> bsize(neighbors.size());
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
					bsize[neigh] += lambdas;
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

	for (size_t i = 0; i < neighbors.size(); ++i) {
		send[i].resize(bsize[i]);
		recv[i].resize(bsize[i]);
	}
}

}