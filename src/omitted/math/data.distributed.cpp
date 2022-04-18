
#include "data.distributed.h"
#include "data.synchronization.h"
#include "data.mv.h"
#include "vector.sparse.h"
#include "esinfo/eslog.h"

#include <cstddef>
#include <cstring>
#include <algorithm>

using namespace espreso;

_DataDistributed::_DataDistributed()
: nhalo(0), ranks(0), nneighbors(0), halo(NULL), distribution(NULL), nintervals(NULL), neighbors(NULL), sync(NULL), mv(NULL)
{

}

void _DataDistributed::alloc(esint nhalo, esint ranks, esint nneighbors)
{
	clear();
	this->nhalo = nhalo;
	this->ranks = ranks;
	this->nneighbors = nneighbors;

	if (nhalo) {
		halo = new esint[nhalo];
	}
	distribution = new esint[ranks + 1];
	nintervals = new esint[2 * nneighbors];
	neighbors = new int[nneighbors];
	sync = new DataSynchronization();
	mv = new DataMV();
}

void _DataDistributed::clear()
{
	nhalo = 0;
	ranks = 0;
	nneighbors = 0;
	if (halo) { delete[] halo; halo = NULL; }
	if (distribution) { delete[] distribution; distribution = NULL; }
	if (nintervals) { delete[] nintervals; nintervals = NULL; }
	if (neighbors) { delete[] neighbors; neighbors = NULL; }
	if (sync) { delete sync; sync = NULL; }
	if (mv) { delete mv; mv = NULL; }
}


DataDistributed::DataDistributed()
{

}

DataDistributed::DataDistributed(esint nhalo, esint ranks, esint nneighbors)
{
	resize(nhalo, ranks, nneighbors);
}

DataDistributed::~DataDistributed()
{
	_allocated.clear();
}

void DataDistributed::resize(esint nhalo, esint ranks, esint nneighbors)
{
	_allocated.clear();
	_allocated.alloc(nhalo, ranks, nneighbors);
	_DataDistributed::operator=(_allocated);
}

void DataDistributed::swap(DataDistributed *other)
{
	_DataDistributed _data = _allocated;
	_allocated = other->_allocated;
	other->_allocated = _data;

	_data = *this;
	this->_DataDistributed::operator=(*other);
	other->_DataDistributed::operator=(_data);
}

void DataDistributed::shallowCopy(const DataDistributed *other)
{
	_allocated.clear();
	_DataDistributed::operator=(*other);
}


void DataDistributed::deepCopy(const DataDistributed *other)
{
	_allocated.clear();
	_allocated.alloc(other->nhalo, other->ranks, other->nneighbors);
	_DataDistributed::operator=(_allocated);
	fillDistribution(other->halo, other->distribution, other->neighbors);
	*sync = *other->sync;
	*mv = *other->mv;
}

void DataDistributed::uniformCombination(const DataDistributed *first, const DataDistributed *second, esint nfirst, esint nsecond)
{
	resize(first->nhalo + second->nhalo, first->ranks, first->nneighbors);
	VectorSparse::combineIndices(halo, first->halo, second->halo, first->halo + first->nhalo, second->halo + second->nhalo, nfirst, nsecond);
	for (int r = 0; r <= ranks; r++) {
		distribution[r] = (nfirst + nsecond) * (first->distribution[r] / nfirst);
	}
	for (int n = 0; n < nneighbors; n++) {
		nintervals[2 * n + 0] = (nfirst + nsecond) * (first->nintervals[2 * n + 0] / nfirst);
		nintervals[2 * n + 1] = (nfirst + nsecond) * (first->nintervals[2 * n + 1] / nfirst);
	}
	memcpy(neighbors, first->neighbors, sizeof(int) * nneighbors);

	sync->uniformCombination(first->sync, second->sync, nfirst, nsecond);
	mv->uniformCombination(first->mv, second->mv, nfirst, nsecond);
}


void DataDistributed::fillDistribution(esint *halo, esint *distribution, int *neighbors)
{
	memcpy(this->halo, halo, sizeof(esint) * nhalo);
	memcpy(this->distribution, distribution, sizeof(esint) * (ranks + 1));
	memcpy(this->neighbors, neighbors, sizeof(int) * nneighbors);

	for (esint n = 0; n < nneighbors; n++) {
		nintervals[2 * n + 0] = std::lower_bound(halo, halo + nhalo, distribution[neighbors[n]]) - halo;
		nintervals[2 * n + 1] = std::lower_bound(halo, halo + nhalo, distribution[neighbors[n] + 1]) - halo;
	}

	sync->neighbors.assign(neighbors, neighbors + nneighbors);
}

