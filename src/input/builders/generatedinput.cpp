
#include "generatedinput.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.h"

#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

GeneratedInput::GeneratedInput(MeshBuilder &meshData, bool needSynchronization)
: Input(meshData)
{
	eslog::startln("BUILDER: BUILD GENERATED MESH", "BUILDER");
	profiler::syncstart("generated_input");

	removeDanglingNodes();
	profiler::synccheckpoint("remove_dangling");
	eslog::checkpointln("BUILDER: DANGLING NODES REMOVED");

	fillNeighbors();
	profiler::synccheckpoint("fill_neighbors");
	eslog::checkpointln("BUILDER: NEIGHBOURS FILLED");

	if (needSynchronization) {
		synchronizeGlobalIndices();
		profiler::synccheckpoint("synchronize_indices");
		eslog::checkpointln("BUILDER: NODES INDICES SYNCHRONIZED");

		sortNodes(true);
		profiler::synccheckpoint("sort_nodes");
		eslog::checkpointln("BUILDER: NODES SORTED");
	}

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
	profiler::syncend("generated_input");
	eslog::endln("BUILDER: REGIONS FILLED");
}

struct __Point__ {

	static constexpr size_t digits = 1000;

	__Point__(): x(0), y(0), z(0), id(-1) {};
	__Point__(const Point &p, esint id): x(p.x), y(p.y), z(p.z), id(id) {};

	bool operator<(const __Point__ &p) const
	{
		if (std::trunc(z * digits) == std::trunc(p.z * digits)) {
			if (std::trunc(y * digits) == std::trunc(p.y * digits)) {
				if (std::trunc(x * digits) == std::trunc(p.x * digits)) {
					return false;
				}
				return x < p.x;
			}
			return y < p.y;
		}
		return z < p.z;
	}

	bool operator==(const __Point__ &p) const
	{
		return !(*this < p) && !(p < *this);
	}

	double x, y, z;
	esint id;
};

void GeneratedInput::removeDanglingNodes()
{
	profiler::syncstart("remove_dangling_nodes");
	std::vector<esint> usedNodes = _meshData.enodes;
	utils::sortAndRemoveDuplicates(usedNodes);

	std::vector<Point> coordinates;
	std::vector<esint> nIDs, ndist, noffset(_meshData.nIDs.size());
	std::vector<int> nranks;

	ndist.push_back(0);
	for (size_t i = 0; i < usedNodes.size(); i++) {
		noffset[usedNodes[i]] = i;
		coordinates.push_back(_meshData.coordinates[usedNodes[i]]);
		nIDs.push_back(_meshData.nIDs[usedNodes[i]]);
		if (_meshData._nranks.size()) {
			nranks.insert(nranks.end(), _meshData._nranks.begin() + _meshData._nrankdist[usedNodes[i]], _meshData._nranks.begin() + _meshData._nrankdist[usedNodes[i] + 1]);
			ndist.push_back(nranks.size());
		}
	}

	for (size_t i = 0; i < _meshData.enodes.size(); i++) {
		_meshData.enodes[i] = noffset[_meshData.enodes[i]];
	}

	for (auto it = _meshData.nregions.begin(); it != _meshData.nregions.end(); ++it) {
		for (size_t i = 0; i < it->second.size(); i++) {
			it->second[i] = noffset[it->second[i]];
		}
	}

	_meshData.nIDs.swap(nIDs);
	_meshData.coordinates.swap(coordinates);
	if (_meshData._nranks.size()) {
		_meshData._nrankdist.swap(ndist);
		_meshData._nranks.swap(nranks);
	}
	profiler::syncend("remove_dangling_nodes");
}


void GeneratedInput::synchronizeGlobalIndices()
{
	auto n2i = [ & ] (size_t neighbor) {
		return std::lower_bound(info::mesh->neighbors.begin(), info::mesh->neighbors.end(), neighbor) - info::mesh->neighbors.begin();
	};

	std::vector<std::vector<__Point__> > sBuffer(info::mesh->neighbors.size());
	std::vector<std::vector<__Point__> > rBuffer(info::mesh->neighbors.size());

	for (size_t n = 0; n < _meshData.nIDs.size(); ++n) {
		if (_meshData._nrankdist[n + 1] - _meshData._nrankdist[n] > 1) {
			if (_meshData._nranks[_meshData._nrankdist[n]] == info::mpi::rank) {
				for (esint r = _meshData._nrankdist[n] + 1; r < _meshData._nrankdist[n + 1]; ++r) {
					sBuffer[n2i(_meshData._nranks[r])].push_back(__Point__(_meshData.coordinates[n], _meshData.nIDs[n]));
				}
			} else {
				sBuffer[n2i(_meshData._nranks[_meshData._nrankdist[n]])].push_back(__Point__(_meshData.coordinates[n], n));
			}
		}
	}

	for (size_t n = 0; n < sBuffer.size(); n++) {
		std::sort(sBuffer[n].begin(), sBuffer[n].end());
		rBuffer[n].resize(sBuffer[n].size());
	}

	if (!Communication::receiveLowerKnownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
		eslog::error("problem while synchronization of global indices.\n");
	}

	for (size_t n = 0; n < info::mesh->neighbors.size(); n++) {
		if (info::mesh->neighbors[n] < info::mpi::rank) {
			for (size_t p = 0; p < rBuffer[n].size(); p++) {
				auto it = std::lower_bound(sBuffer[n].begin(), sBuffer[n].end(), rBuffer[n][p]);
				if (*it == rBuffer[n][p]) {
					_meshData.nIDs[it->id] = rBuffer[n][p].id;
				} else {
					eslog::error("Internal ERROR while synchronization global indices.\n");
				}
			}
		}
	}
}


