
#include "distributedinfo.h"

#include "../basis/utilities/utils.h"
#include "../configuration/environment.h"

#include "../basis/point/point.h"
#include "../mesh/elements/element.h"
#include "../mesh/settings/property.h"
#include "../mesh/structures/coordinates.h"
#include "../mesh/structures/mesh.h"
#include "../mesh/structures/region.h"

#include <numeric>

using namespace espreso::output;

static std::vector<espreso::Point> computeDomainsCenters(const espreso::Mesh *mesh)
{
	std::vector<espreso::Point> domainsCenters(mesh->parts());
	for (size_t p = 0; p < mesh->coordinates().parts(); p++) {
		for (size_t i = 0; i < mesh->coordinates().localSize(p); i++) {
			domainsCenters[p] += mesh->coordinates().get(i, p);
		}
		domainsCenters[p] /= mesh->coordinates().localSize(p);
	}
	return domainsCenters;
}

static espreso::Point computeClusterCenter(const espreso::Mesh *mesh)
{
	espreso::Point clusterCenter;
	for (size_t i = 0; i < mesh->coordinates().clusterSize(); i++) {
		clusterCenter += mesh->coordinates()[i];
	}
	clusterCenter /= mesh->coordinates().clusterSize();
	return clusterCenter;
}

DistributedInfo::DistributedInfo(const Mesh *mesh, double domainShrinkRatio, double clusterShrinkRatio, InfoMode mode)
: MeshInfo(mesh, mode), _domainShrinkRatio(domainShrinkRatio), _clusterShrinkRatio(clusterShrinkRatio)
{
	if (_mode | InfoMode::PREPARE) {
		_domainsCenters = computeDomainsCenters(_mesh);
		_clusterCenter = computeClusterCenter(_mesh);
		prepare(_mesh->elements(), 0, _mesh->elements().size());
	}
}

DistributedInfo::DistributedInfo(const Mesh *mesh, const Region* region, double domainShrinkRatio, double clusterShrinkRatio, InfoMode mode)
: MeshInfo(mesh, region, mode), _domainShrinkRatio(domainShrinkRatio), _clusterShrinkRatio(clusterShrinkRatio)
{
	_domainsCenters = computeDomainsCenters(_mesh);
	_clusterCenter = computeClusterCenter(_mesh);
	prepare(region->elements(), 0, region->elements().size());
}


MeshInfo* DistributedInfo::deriveRegion(const Region *region) const
{
	DistributedInfo *copy = new DistributedInfo(_mesh, _domainShrinkRatio, _clusterShrinkRatio);
	copy->_domainsCenters = _domainsCenters;
	copy->_clusterCenter = _clusterCenter;
	copy->_region = region;
	copy->prepare(region->elements(), 0, region->elements().size());
	return copy;
}

MeshInfo* DistributedInfo::copyWithoutMesh() const
{
	DistributedInfo *copy = new DistributedInfo(_mesh, _domainShrinkRatio, _clusterShrinkRatio);
	copy->_domainsCenters = _domainsCenters;
	copy->_clusterCenter = _clusterCenter;
	copy->_regions.push_back(RegionData());
	return copy;
}

espreso::Point DistributedInfo::shrink(const Point &p, eslocal domain) const
{
	Point point = _clusterCenter + (p - _clusterCenter) * _clusterShrinkRatio;
	point = _domainsCenters[domain] + (point - _domainsCenters[domain]) * _domainShrinkRatio;
	return point;
}

void DistributedInfo::prepare(const std::vector<Element*> &region, size_t begin, size_t end)
{
	_regions.push_back(RegionData());
	if (!region.size()) {
		return;
	}

	std::vector<std::vector<eslocal> > rCoordinates(_mesh->parts());
	for (size_t e = 0; e < region.size(); e++) {
		for (auto d = region[e]->domains().begin(); d != region[e]->domains().end(); ++d) {
			for (size_t n = 0; n < region[e]->nodes(); n++) {
				rCoordinates[*d].push_back(region[e]->node(n));
			}
		}
	}

	#pragma omp parallel for
	for (size_t p = 0; p < _mesh->parts(); p++) {
		std::sort(rCoordinates[p].begin(), rCoordinates[p].end());
		Esutils::removeDuplicity(rCoordinates[p]);
	}

	std::vector<size_t> offsets;
	// compute coordinates
	std::vector<std::vector<double> > dCoordinates(_mesh->parts());
	#pragma omp parallel for
	for (size_t p = 0; p < rCoordinates.size(); p++) {
		Point point;
		dCoordinates[p].reserve(rCoordinates[p].size());
		for (size_t i = 0; i < rCoordinates[p].size(); i++) {
			point = _mesh->coordinates()[rCoordinates[p][i]];
			point = _clusterCenter + (point - _clusterCenter) * _clusterShrinkRatio;
			point = _domainsCenters[p] + (point - _domainsCenters[p]) * _domainShrinkRatio;
			dCoordinates[p].insert(dCoordinates[p].end(), { point.x, point.y, point.z });
		}
	}

	_regions.back().coordinates.clear();
	offsets = { 0 };
	for (size_t p = 0; p < _mesh->parts(); p++) {
		_regions.back().coordinates.insert(_regions.back().coordinates.end(), dCoordinates[p].begin(), dCoordinates[p].end());
		offsets.push_back(_regions.back().coordinates.size() / 3);
	}

	for (size_t e = 0, offset = 0; e < region.size(); e++) {
		_regions.back().elementsTypes.insert(_regions.back().elementsTypes.end(), region[e]->domains().size(), region[e]->vtkCode());
		for (auto d = region[e]->domains().begin(); d != region[e]->domains().end(); ++d, offset += region[e]->nodes()) {
			_regions.back().elementsNodes.push_back(offset + region[e]->nodes());
			for (size_t n = 0; n < region[e]->nodes(); n++) {
				eslocal oIndex = std::lower_bound(rCoordinates[*d].begin(), rCoordinates[*d].end(), region[e]->node(n)) - rCoordinates[*d].begin();
				_regions.back().elements.push_back(oIndex + offsets[*d]);
			}
		}
	}
}

void DistributedInfo::addGeneralInfo()
{
	if (_region == NULL) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: general info can be added only to region with all elements.";
	}

	std::vector<eslocal> *pointIDcluster = new std::vector<eslocal>(_regions.back().coordinates.size() / 3);
	std::vector<eslocal> *pointIDglobal = new std::vector<eslocal>();
	std::vector<eslocal> *elementID = new std::vector<eslocal>(_mesh->getPartition().back());
	std::vector<eslocal> *decomposition = new std::vector<eslocal>();

	std::iota(pointIDcluster->begin(), pointIDcluster->end(), 0);
	std::iota(elementID->begin(), elementID->end(), 0);
	for (size_t p = 0; p < _mesh->getPartition().size() - 1; p++) {
		decomposition->insert(decomposition->end(), _mesh->getPartition()[p + 1] - _mesh->getPartition()[p], p);
	}

	pointIDglobal->reserve(_regions.back().coordinates.size() / 3);
	for (size_t p = 0; p < _mesh->parts(); p++) {
		for (size_t i = 0; i < _mesh->coordinates().localSize(p); i++) {
			pointIDglobal->push_back(_mesh->coordinates().globalIndex(i, p));
		}
	}

	_regions.back().data.pointDataInteger["pointIDcluster"] = std::make_pair(1, pointIDcluster);
	_regions.back().data.pointDataInteger["pointIDglobal"] = std::make_pair(1, pointIDglobal);
	_regions.back().data.elementDataInteger["elementID"] = std::make_pair(1, elementID);
	_regions.back().data.elementDataInteger["decomposition"] = std::make_pair(1, decomposition);

	for (int r = 0; r < espreso::environment->MPIsize; r++) {
		std::vector<eslocal> *cluster = new std::vector<eslocal>(_regions.back().coordinates.size() / 3);
		_regions.back().data.pointDataInteger["cluster" + std::to_string(r)] = std::make_pair(1, cluster);
	}

	std::vector<size_t> offsets(_mesh->parts());
	for (size_t p = 1; p < _mesh->parts(); p++) {
		offsets[p] = offsets[p - 1] + _mesh->coordinates().localSize(p - 1);
	}

	for (size_t n = 0; n < _mesh->nodes().size(); n++) {
		for (auto c = _mesh->nodes()[n]->clusters().begin(); c != _mesh->nodes()[n]->clusters().end(); ++c) {
			for (auto d = _mesh->nodes()[n]->domains().begin(); d != _mesh->nodes()[n]->domains().end(); ++d) {
				(*_regions.back().data.pointDataInteger["cluster" + std::to_string(*c)].second)[offsets[*d] + _mesh->coordinates().localIndex(n, *d)] = 1;
			}
		}
	}
}

void DistributedInfo::addSettings(size_t step)
{
	const Region *region = _region;
	if (region == NULL) {
		region = _mesh->regions()[0];
	}

	if (region->settings.size() <= step) {
		return;
	}
	for (auto it = region->settings[step].begin(); it != region->settings[step].end(); ++it) {
		std::vector<double> *values = new std::vector<double>();

		values->reserve(region->elements().size());
		for (size_t e = 0; e < region->elements().size(); e++) {
			values->push_back(0);
			for (size_t n = 0; n < region->elements()[e]->nodes(); n++) {
				values->back() += region->elements()[e]->sumProperty(it->first, n, step, 0);
			}
			values->back() /= region->elements()[e]->nodes();
			values->insert(values->end(), region->elements()[e]->domains().size() - 1, values->back());
		}

		std::stringstream ss;
		ss << it->first;
		_regions.back().data.elementDataDouble[ss.str()] = std::make_pair(1, values);
	}
}

void DistributedInfo::addSolution(const std::vector<Solution*> &solution)
{
	_regions.back().solutions.insert(_regions.back().solutions.end(), solution.begin(), solution.end());
}


