
#include "sequentialinput.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"

#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

SequentialInput::SequentialInput(MeshBuilder &mesh)
: Input(mesh)
{
	eslog::startln("BUILDER: BUILD SEQUENTIAL MESH", "BUILDER");
	profiler::syncstart("sequential_input");

	if (_meshData.removeDuplicates) {
		searchDuplicateNodes();
		coupleDuplicateNodes();
		profiler::synccheckpoint("merge_duplicate_nodes");
	}

	sortNodes(true);
	removeDanglingNodes();
	_nDistribution = { 0, (esint)_meshData.nIDs.size() };
	profiler::synccheckpoint("sort_nodes");
	eslog::checkpointln("BUILDER: NODES SORTED");

	_eDistribution = { 0, (esint)_meshData.eIDs.size() };
	sortElements();
	profiler::synccheckpoint("sort_elements");
	eslog::checkpointln("BUILDER: ELEMENTS SORTED");

	if (_meshData.removeDuplicates) {
		removeDuplicateElements();
		profiler::synccheckpoint("merge_duplicate_elements");
	}

	if (_meshData._duplicateElements.size()) {
		for (auto it = _meshData.eregions.begin(); it != _meshData.eregions.end(); ++it) {
			std::vector<esint> rdata;
			auto remove = _meshData._duplicateElements.begin();
			for (auto e = it->second.begin(); e != it->second.end(); ++e) {
				while (remove != _meshData._duplicateElements.end() && remove->id < *e) { ++remove; }
				if (remove == _meshData._duplicateElements.end() || remove->id != *e) {
					rdata.push_back(*e);
				}
				if (remove != _meshData._duplicateElements.end() && remove->id == *e) {
					rdata.push_back(remove->target);
				}
			}
			utils::sortAndRemoveDuplicates(rdata);
			it->second = rdata;
		}
		profiler::synccheckpoint("merge_duplicate_eregions");
	}

	if (_meshData._duplicateNodes.size()) {
		for (auto it = _meshData.nregions.begin(); it != _meshData.nregions.end(); ++it) {
			std::vector<esint> rdata;
			auto projection = _meshData._duplicateNodes.begin();
			for (auto n = it->second.begin(), nit = _meshData.nIDs.begin(); n != it->second.end(); ++n) {
				while (nit != _meshData.nIDs.end() && *nit < *n) { ++nit; }
				if (nit != _meshData.nIDs.end() && *nit == *n) {
					rdata.push_back(*n);
				}
				while (projection != _meshData._duplicateNodes.end() && projection->id < *n) { ++projection; }
				if (projection != _meshData._duplicateNodes.end() && projection->id == *n) {
					rdata.push_back(projection->target);
				}
			}
			utils::sortAndRemoveDuplicates(rdata);
			it->second = rdata;
		}
		profiler::synccheckpoint("merge_duplicate_nregions");
	}

	reindexERegions();
	reindexNRegions();
	profiler::synccheckpoint("reindex_regions");
	eslog::checkpointln("BUILDER: REGION REINDEXED");

	_meshData._nranks.resize(_meshData.nIDs.size());
	_meshData._nrankdist.resize(_meshData.nIDs.size() + 1);
	std::iota(_meshData._nrankdist.begin(), _meshData._nrankdist.end(), 0);
	info::mesh->neighborsWithMe.push_back(0);
	profiler::synccheckpoint("fill_rank");
	eslog::checkpointln("BUILDER: RANKS FILLED");

	clip();

	fillNodes();
	profiler::synccheckpoint("fill_nodes");
	eslog::checkpointln("BUILDER: NODES FILLED");

	fillElements();
	profiler::synccheckpoint("fill_elements");
	eslog::checkpointln("BUILDER: ELEMENTS FILLED");

	fillElementRegions();
	fillBoundaryRegions();
	fillNodeRegions();
	profiler::synccheckpoint("fill_regions");
	eslog::checkpointln("BUILDER: REGIONS FILLED");
	profiler::syncend("sequential_input");
	eslog::endln("BUILDER: ELEMENTS NODES REINDEXED");
}

void SequentialInput::removeDanglingNodes()
{
	profiler::syncstart("remove_dangling_nodes");
	std::vector<esint> usedNodes; usedNodes.reserve(_meshData.enodes.size());
	for (size_t e = 0, offset = 0; e < _meshData.esize.size(); offset += _meshData.esize[e++]) {
		if (_meshData.etype[e] != (int)Element::CODE::NOT_SUPPORTED) {
			usedNodes.insert(usedNodes.end(), _meshData.enodes.begin() + offset, _meshData.enodes.begin() + offset + _meshData.esize[e]);
		}
	}
	utils::sortAndRemoveDuplicates(usedNodes);

	std::vector<Point> coordinates;
	std::vector<esint> nIDs, ndist, noffset(usedNodes.back() + 1);
	std::vector<int> nranks;
	coordinates.reserve(_meshData.nIDs.size());
	nIDs.reserve(_meshData.nIDs.size());
	ndist.reserve(_meshData.nIDs.size() + 1);
	nranks.reserve(_meshData._nranks.size());

	ndist.push_back(0);
	for (size_t n = 0, i = 0; n < _meshData.nIDs.size(); ++n) {
		if (i < usedNodes.size() && usedNodes[i] == _meshData.nIDs[n]) {
			noffset[usedNodes[i]] = i;
			coordinates.push_back(_meshData.coordinates[n]);
			nIDs.push_back(_meshData.nIDs[n]);
			if (_meshData._nranks.size()) {
				nranks.insert(nranks.end(), _meshData._nranks.begin() + _meshData._nrankdist[n], _meshData._nranks.begin() + _meshData._nrankdist[n + 1]);
				ndist.push_back(nranks.size());
			}
			++i;
		} else {
			_dangling.push_back(_meshData.nIDs[n]);
		}
	}

	for (size_t i = 0; i < _meshData.enodes.size(); i++) {
		_meshData.enodes[i] = noffset[_meshData.enodes[i]];
	}

	_meshData.nIDs.swap(nIDs);
	_meshData.coordinates.swap(coordinates);
	if (_meshData._nranks.size()) {
		_meshData._nrankdist.swap(ndist);
		_meshData._nranks.swap(nranks);
	}
	profiler::syncend("remove_dangling_nodes");
}

void SequentialInput::reindexNRegions()
{
	profiler::syncstart("reindex_nregions");
	for (auto region = _meshData.nregions.begin(); region != _meshData.nregions.end(); ++region) {
		auto n = region->second.begin();
		if (n != region->second.end()) {
			auto nit = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.begin(), *n);
			size_t i = 0;
			for (; n != region->second.end(); ++n) {
				while (nit != _meshData.nIDs.end() && *nit < *n) { ++nit; }
				if (nit == _meshData.nIDs.end() || *n != *nit) {
					if (!std::binary_search(_dangling.begin(), _dangling.end(), *n)) {
						eslog::error("MESIO: node region '%s' contains ID[%d] that is not found\n", region->first.c_str(), *n);
					}
				} else {
					region->second[i++] = nit - _meshData.nIDs.begin();
				}
			}
			region->second.resize(i);
		}
	}
	profiler::syncend("reindex_nregions");
}

void SequentialInput::reindexERegions()
{
	profiler::syncstart("reindex_eregions");
	for (auto region = _meshData.eregions.begin(); region != _meshData.eregions.end(); ++region) {
		auto n = region->second.begin();
		if (n != region->second.end()) {
			auto nit = _meshData.eIDs.begin();
			for (size_t t = 0, begin = 0; t < _etypeDistribution.size(); ++t) {
				if (nit != _meshData.eIDs.end() && *nit != *n) {
					nit = std::lower_bound(_meshData.eIDs.begin() + begin, _meshData.eIDs.begin() + _etypeDistribution[t], *n);
					begin = _etypeDistribution[t];
				}
			}
			for (; n != region->second.end(); ++n) {
				while (nit != _meshData.eIDs.end() && *nit < *n) { ++nit; }
				if (nit == _meshData.eIDs.end() || *n != *nit) {
					eslog::error("MESIO: element region '%s' contains ID[%d] that is not found\n", region->first.c_str(), *n);
				} else {
					*n = nit - _meshData.eIDs.begin();
				}
			}
		}
	}
	profiler::syncend("reindex_eregions");
}



