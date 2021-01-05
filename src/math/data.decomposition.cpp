
#include "data.decomposition.h"
#include "vector.dense.feti.h"
#include "esinfo/eslog.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "basis/containers/serializededata.h"
#include "wrappers/mpi/communication.h"

#include <cstddef>
#include <vector>
#include <cstring>

using namespace espreso;

struct espreso::DataDecompositionGather {
	std::vector<std::vector<esint> > sindices, rindices;
	std::vector<std::vector<double> > send, recv;
	std::vector<int> neighbors;
	esint gatheredSize;
};

_DataDecomposition::_DataDecomposition()
: rank(0), ranks(1), nneighbors(0), distribution(NULL), neighbors(NULL), dmap(NULL),
  doffset(0), nshared(0), shared(NULL),
  gathered(NULL), gather(NULL)
{

}

void _DataDecomposition::alloc(esint ranks, esint nneighbors)
{
	clear();
	this->ranks = ranks;
	this->nneighbors = nneighbors;

	distribution = new esint[ranks + 1];
	neighbors = new int[nneighbors];
	gather = new DataDecompositionGather();
}

void _DataDecomposition::clear()
{
	rank = 0;
	ranks = 1;
	nneighbors = 0;
	if (distribution) { delete[] distribution; distribution = NULL; }
	if (neighbors) { delete[] neighbors; neighbors = NULL; }
	if (dmap) { delete dmap; dmap = NULL; }

	doffset = 0;
	nshared = 0;
	if (shared) { delete[] shared; shared = NULL; }

	if (gathered) { delete[] gathered; gathered = NULL; }
	if (gather) { delete gather; gather = NULL; }
}

DataDecomposition::DataDecomposition(DUPLICATION duplications)
: duplications(duplications)
{

}

DataDecomposition::~DataDecomposition()
{
	_allocated.clear();
}

void DataDecomposition::swap(DataDecomposition *other)
{
	_DataDecomposition _data = _allocated;
	_allocated = other->_allocated;
	other->_allocated = _data;

	_data = *this;
	this->_DataDecomposition::operator=(*other);
	other->_DataDecomposition::operator=(_data);

	DUPLICATION d = duplications;
	duplications = other->duplications;
	other->duplications = d;
}

void DataDecomposition::shallowCopy(const DataDecomposition *other)
{
	_allocated.clear();
	_DataDecomposition::operator=(*other);
	duplications = other->duplications;
}

void DataDecomposition::shallowCopyStructure(const DataDecomposition *other)
{
	shallowCopy(other);
	if (gather != NULL) {
		_allocated.gathered = gathered = new double[gather->gatheredSize];
	}
}

void DataDecomposition::deepCopy(const DataDecomposition *other)
{
	deepCopyStructure(other);
	if (gather != NULL) {
		memcpy(gathered, other->gathered, sizeof(double) * gather->gatheredSize);
	}
}

void DataDecomposition::deepCopyStructure(const DataDecomposition *other)
{
	_allocated.clear();
	_allocated.alloc(other->ranks, other->nneighbors);
	_DataDecomposition::operator=(_allocated);
	fillDecomposition(other->rank, other->ranks, other->nneighbors, other->distribution, other->neighbors, other->dmap);
	duplications = other->duplications;
	_allocated.gather = new DataDecompositionGather(*other->gather);
	_allocated.gathered = gathered = new double[gather->gatheredSize];
}

void DataDecomposition::uniformCombination(const DataDecomposition *first, const DataDecomposition *second, int nfirst, int nsecond)
{
	_allocated.clear();
	_allocated.alloc(first->ranks, first->nneighbors);
	_DataDecomposition::operator=(_allocated);
	this->rank = first->rank;
	duplications = first->duplications;
	memcpy(this->distribution, first->distribution, sizeof(esint) * (ranks + 1));
	memcpy(this->neighbors, first->neighbors, sizeof(int) * nneighbors);
	doffset = distribution[rank];
	_allocated.dmap = dmap = combineDomainMap(first->dmap, second->dmap, nfirst, nsecond);

	buildGatherData();
}

void DataDecomposition::fillDecomposition(esint rank, esint ranks, esint nneighbors, esint *distribution, int *neighbors, const serializededata<esint, DI> *dmap)
{
	_allocated.clear();
	_allocated.alloc(ranks, nneighbors);
	_DataDecomposition::operator=(_allocated);
	this->rank = rank;
	memcpy(this->distribution, distribution, sizeof(esint) * (ranks + 1));
	memcpy(this->neighbors, neighbors, sizeof(int) * nneighbors);
	doffset = distribution[rank];
	_allocated.dmap = this->dmap = new serializededata<esint, DI>(*dmap);

	buildGatherData();
}

void DataDecomposition::allGather(const VectorDenseFETI &in) const
{
	auto map = dmap->begin();
	for (esint n = 0, i = 0, prev = 0; n < nshared; prev = shared[n++]) {
		map += shared[n] - prev;
		for (auto di = map->begin(); di != map->end(); ++di, ++i) {
			if (ismy(di->domain)) {
				gathered[i] = in[di->domain - doffset][di->index];
			}
		}
	}

	for (size_t n = 0; n < gather->sindices.size(); ++n) {
		for (size_t i = 0; i < gather->sindices[n].size(); ++i) {
			gather->send[n][i] = gathered[gather->sindices[n][i]];
		}
	}

	if (!Communication::exchangeKnownSize(gather->send, gather->recv, gather->neighbors)) {
		eslog::internalFailure("gather decomposed data.\n");
	}

	for (size_t n = 0; n < gather->rindices.size(); ++n) {
		for (size_t i = 0; i < gather->rindices[n].size(); ++i) {
			gathered[gather->rindices[n][i]] = gather->recv[n][i];
		}
	}
}

serializededata<esint, DI>* DataDecomposition::combineDomainMap(const serializededata<esint, DI> *first, const serializededata<esint, DI> *second, int nfirst, int nsecond)
{
	std::vector<esint> distribution({ 0 });
	std::vector<DI> data;
	std::vector<size_t> tdistribution = first->boundarytarray().distribution();
	std::vector<size_t> tdata = first->datatarray().distribution();

	int nsum = nfirst + nsecond;
	auto dfirst = first->begin();
	auto dsecond = second->begin();
	for (size_t n = 0; n < first->structures() / nfirst; ++n) {
		for (int i = 0; i < nfirst; ++i, ++dfirst) {
			for (auto di = dfirst->begin(); di != dfirst->end(); ++di) {
				data.push_back({ di->domain, nsum * (di->index / nfirst) + di->index % nfirst });
			}
			distribution.push_back(data.size());
		}
		for (int j = 0; j < nsecond; ++j, ++dsecond) {
			for (auto di = dsecond->begin(); di != dsecond->end(); ++di) {
				data.push_back({ di->domain, nsum * (di->index / nsecond) + di->index % nsecond + nfirst});
			}
			distribution.push_back(data.size());
		}
	}

	for (size_t i = 0; i < tdistribution.size(); ++i) {
		tdistribution[i] += second->boundarytarray().distribution()[i] - (i ? 1 : 0);
		tdata[i] += second->datatarray().distribution()[i];
		if (tdistribution[i] > distribution.size() + 1) {
			tdistribution[i] = distribution.size() + 1;
		}
		if (tdata[i] > data.size() + 1) {
			tdata[i] = data.size() + 1;
		}
	}

	return new serializededata<esint, DI>(
				tarray<esint>(tdistribution, distribution),
				tarray<DI>(tdata, data));
}

void DataDecomposition::buildGatherData()
{
	gather->neighbors.assign(neighbors, neighbors + nneighbors);
	gather->sindices.resize(gather->neighbors.size());
	gather->rindices.resize(gather->neighbors.size());
	gather->send.resize(gather->neighbors.size());
	gather->recv.resize(gather->neighbors.size());

	std::vector<esint> shindices;
	esint ii = 0;
	for (auto map = dmap->begin(); map != dmap->end(); ++map, ++ii) {
		if (map->size() > 1) {
			shindices.push_back(ii);
			gather->gatheredSize += map->size();
		}
	}
	nshared = shindices.size();
	_allocated.shared = shared = new esint[nshared];
	memcpy(shared, shindices.data(), sizeof(esint) * shindices.size());

	if (nshared) {
		esint i = 0, offset = 0, size = 0, prev = 0;
		auto map = dmap->begin();
		for (esint n = 0; n < nshared; i += map->size(), prev = shared[n++]) {
			map += shared[n] - prev;
			size = 0;
			for (auto di = map->begin(); di != map->end(); ++di) {
				if (ismy(di->domain)) {
					if (size++ == 0) {
						offset = di - map->begin();
					}
				}
			}
			auto current = map->begin();
			for (esint n = 0; n < nneighbors; ++n) {
				while (current != map->end() && ismy(current->domain)) {
					++current;
				}
				auto begin = current;
				while (current != map->end() && current->domain < distribution[neighbors[n] + 1]) {
					++current;
				}
				if (begin != map->end() && !ismy(begin->domain)) {
					if (begin != current) {
						for (esint j = i + offset; j < i + offset + size; ++j) {
							gather->sindices[n].push_back(j);
						}
					}
					for (auto di = begin; di != current; ++di) {
						gather->rindices[n].push_back(i + (di - map->begin()));
					}
				}
			}
		}

		for (size_t n = 0; n < gather->neighbors.size(); ++n) {
			gather->send[n].resize(gather->sindices[n].size());
			gather->recv[n].resize(gather->rindices[n].size());
		}
	}

	_allocated.gathered = gathered = new double[gather->gatheredSize];
}

