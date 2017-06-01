
#include "distributedinfo.h"

#include "../basis/utilities/utils.h"
#include "../configuration/environment.h"

#include "../basis/point/point.h"
#include "../mesh/elements/element.h"
#include "../mesh/settings/property.h"
#include "../mesh/structures/coordinates.h"
#include "../mesh/structures/mesh.h"
#include "../mesh/structures/elementtypes.h"
#include "../mesh/structures/region.h"

#include "../assembler/solution.h"

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
	if (_mode & InfoMode::PREPARE) {
		_domainsCenters = computeDomainsCenters(_mesh);
		_clusterCenter = computeClusterCenter(_mesh);
		prepare(_mesh->elements(), _mode);
	}
}

DistributedInfo::DistributedInfo(const Mesh *mesh, const Region* region, double domainShrinkRatio, double clusterShrinkRatio, InfoMode mode)
: MeshInfo(mesh, region, mode), _domainShrinkRatio(domainShrinkRatio), _clusterShrinkRatio(clusterShrinkRatio)
{
	if (_mode & InfoMode::PREPARE) {
		_domainsCenters = computeDomainsCenters(_mesh);
		_clusterCenter = computeClusterCenter(_mesh);
		prepare(region->elements(), _mode);
	}
}


MeshInfo* DistributedInfo::deriveRegion(const Region *region) const
{
	DistributedInfo *copy = new DistributedInfo(_mesh, _domainShrinkRatio, _clusterShrinkRatio, _mode & ~MeshInfo::PREPARE);
	copy->_domainsCenters = _domainsCenters;
	copy->_clusterCenter = _clusterCenter;
	copy->_region = region;
	copy->prepare(region->elements(), copy->_mode);
	return copy;
}

MeshInfo* DistributedInfo::copyWithoutMesh() const
{
	DistributedInfo *copy = new DistributedInfo(_mesh, _domainShrinkRatio, _clusterShrinkRatio, _mode & ~MeshInfo::PREPARE);
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

void DistributedInfo::prepare(const std::vector<Element*> &region, InfoMode mode)
{
	size_t bodies    = mode & InfoMode::SEPARATE_BODIES    ? _mesh->bodies()           : 1;
	size_t materials = mode & InfoMode::SEPARATE_MATERIALS ? _mesh->materials().size() : 1;

	_regions.resize(bodies * materials);
	if (!region.size()) {
		return;
	}

	_cIndices = std::vector<std::vector<std::vector<eslocal> > >(_regions.size(), std::vector<std::vector<eslocal> >(_mesh->parts()));

	size_t threads = environment->OMP_NUM_THREADS;

	// (body * material) x domain x threads x indices
	std::vector<std::vector<std::vector<std::vector<eslocal> > > >
	cIndices(bodies * materials, std::vector<std::vector<std::vector<eslocal> > >(_mesh->parts(), std::vector<std::vector<eslocal> >(threads))),
	eIndices(bodies * materials, std::vector<std::vector<std::vector<eslocal> > >(_mesh->parts(), std::vector<std::vector<eslocal> >(threads))),
	eNodes  (bodies * materials, std::vector<std::vector<std::vector<eslocal> > >(_mesh->parts(), std::vector<std::vector<eslocal> >(threads))),
	eTypes  (bodies * materials, std::vector<std::vector<std::vector<eslocal> > >(_mesh->parts(), std::vector<std::vector<eslocal> >(threads)));

	std::vector<size_t> distribution = Esutils::getDistribution(threads, region.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			size_t regionOffset = 0;
			if ((mode & InfoMode::SEPARATE_BODIES) && region[e]->params()) {
				regionOffset += region[e]->param(Element::Params::BODY) * materials;
			}
			if ((mode & InfoMode::SEPARATE_MATERIALS) && region[e]->params()) {
				regionOffset += region[e]->param(Element::Params::MATERIAL);
			}

			for (auto d = region[e]->domains().begin(); d != region[e]->domains().end(); ++d) {
				eNodes[regionOffset][*d][t].push_back(region[e]->nodes());
				eTypes[regionOffset][*d][t].push_back(region[e]->vtkCode());
				for (size_t n = 0; n < region[e]->nodes(); n++) {
					cIndices[regionOffset][*d][t].push_back(region[e]->node(n));
					eIndices[regionOffset][*d][t].push_back(region[e]->node(n));
				}
			}
		}
	}

	#pragma omp parallel for
	for (size_t p = 0; p < _mesh->parts(); p++) {
		for (size_t r = 0; r < _regions.size(); r++) {
			for (size_t t = 1; t < threads; t++) {
				cIndices[r][p][0].insert(cIndices[r][p][0].end(), cIndices[r][p][t].begin(), cIndices[r][p][t].end());
			}
			std::sort(cIndices[r][p][0].begin(), cIndices[r][p][0].end());
			Esutils::removeDuplicity(cIndices[r][p][0]);
			_cIndices[r][p] = cIndices[r][p][0];
		}
	}

	// (body * material) x domain x coordinates
	std::vector<std::vector<std::vector<double> > > cValues(
			bodies * materials, std::vector<std::vector<double> >(_mesh->parts())
	);

	_cOffset.resize(_regions.size());
	for (size_t r = 0; r < _regions.size(); r++) {

		#pragma omp parallel for
		for (size_t p = 0; p < _mesh->parts(); p++) {
			for (size_t i = 0; i < _cIndices[r][p].size(); i++) {
				Point point = shrink(_mesh->coordinates()[_cIndices[r][p][i]], p);
				cValues[r][p].insert(cValues[r][p].end(), { point.x, point.y, point.z });
			}
		}

		_cOffset[r].resize(_mesh->parts());
		for (size_t p = 0; p < _mesh->parts(); p++) {
			_regions[r].coordinates.insert(_regions[r].coordinates.end(), cValues[r][p].begin(), cValues[r][p].end());
			_cOffset[r][p] = cValues[r][p].size() / 3;
		}
		Esutils::sizesToOffsets(_cOffset[r]);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t p = 0; p < _mesh->parts(); p++) {
				for (size_t i = 0; i < eIndices[r][p][t].size(); i++) {
					eIndices[r][p][t][i] = _cOffset[r][p] + std::lower_bound(_cIndices[r][p].begin(), _cIndices[r][p].end(), eIndices[r][p][t][i]) - _cIndices[r][p].begin();
				}
			}
		}

	}

	#pragma omp parallel for
	for (size_t r = 0; r < _regions.size(); r++) {
		for (size_t p = 0; p < _mesh->parts(); p++) {
			for (size_t t = 0; t < threads; t++) {
				_regions[r].elementsTypes.insert(_regions[r].elementsTypes.end(), eTypes  [r][p][t].begin(), eTypes  [r][p][t].end());
				_regions[r].elementsNodes.insert(_regions[r].elementsNodes.end(), eNodes  [r][p][t].begin(), eNodes  [r][p][t].end());
				_regions[r].elements     .insert(_regions[r].elements.end()     , eIndices[r][p][t].begin(), eIndices[r][p][t].end());
			}
		}
	}


	for (size_t r = 0; r < _regions.size(); r++) {
		distribution = Esutils::getDistribution(threads, _regions[r].elementsNodes.size());

		std::vector<size_t> offsets(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			size_t offset = 0;
			for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
				offset += _regions[r].elementsNodes[e];
			}
			offsets[t] = offset;
		}
		Esutils::sizesToOffsets(offsets);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = distribution[t], offset = offsets[t]; e < distribution[t + 1]; e++) {
				offset = _regions[r].elementsNodes[e] += offset;
			}
		}
	}
}

void DistributedInfo::addGeneralInfo()
{
	if (_region != NULL) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: general info can be added only to region with all elements.";
	}

	size_t materials = _mode & InfoMode::SEPARATE_MATERIALS ? _mesh->materials().size() : 1;

	std::vector<std::vector<eslocal>*> pointIDcluster;
	std::vector<std::vector<eslocal>*> pointIDglobal;
	std::vector<std::vector<eslocal>*> elementID;

	for (size_t r = 0; r < _regions.size(); r++) {
		pointIDcluster.push_back(new std::vector<eslocal>(_regions[r].coordinates.size() / 3));
		pointIDglobal.push_back(new std::vector<eslocal>(_regions[r].coordinates.size() / 3));
		elementID.push_back(new std::vector<eslocal>());

		for (int rank = 0; rank < espreso::environment->MPIsize; rank++) {
			std::vector<eslocal> *cluster = new std::vector<eslocal>(_regions[r].coordinates.size() / 3);
			_regions[r].data.pointDataInteger["cluster" + std::to_string(rank)] = std::make_pair(1, cluster);
		}
	}

	const std::vector<Element*> &region = _region != NULL ? _region->elements() : _mesh->elements();

	for (size_t e = 0; e < region.size(); e++) {
		size_t regionOffset = 0;
		if ((_mode & InfoMode::SEPARATE_BODIES) && region[e]->params()) {
			regionOffset += region[e]->param(Element::Params::BODY) * materials;
		}
		if ((_mode & InfoMode::SEPARATE_MATERIALS) && region[e]->params()) {
			regionOffset += region[e]->param(Element::Params::MATERIAL);
		}

		for (auto d = region[e]->domains().begin(); d != region[e]->domains().end(); ++d) {
			elementID[regionOffset]->push_back(e);
			for (size_t n = 0; n < region[e]->nodes(); n++) {
				eslocal index = std::lower_bound(_cIndices[regionOffset][*d].begin(), _cIndices[regionOffset][*d].end(), region[e]->node(n)) - _cIndices[regionOffset][*d].begin();
				(*pointIDcluster[regionOffset])[index + _cOffset[regionOffset][*d]] = region[e]->node(n);
				(*pointIDglobal[regionOffset])[index + _cOffset[regionOffset][*d]] = _mesh->coordinates().globalIndex((*pointIDcluster[regionOffset])[index + _cOffset[regionOffset][*d]]);

				for (auto c = _mesh->nodes()[region[e]->node(n)]->clusters().begin(); c != _mesh->nodes()[region[e]->node(n)]->clusters().end(); ++c) {
					(*_regions[regionOffset].data.pointDataInteger["cluster" + std::to_string(*c)].second)[index + _cOffset[regionOffset][*d]] = 1;
				}
			}
		}
	}

	for (size_t r = 0; r < _regions.size(); r++) {
		_regions[r].data.pointDataInteger["pointIDcluster"] = std::make_pair(1, pointIDcluster[r]);
		_regions[r].data.pointDataInteger["pointIDglobal"] = std::make_pair(1, pointIDglobal[r]);
		_regions[r].data.elementDataInteger["elementID"] = std::make_pair(1, elementID[r]);
	}
}

void DistributedInfo::addSettings(size_t step)
{
	const Region *region = _region;
	if (region == NULL) {
		region = _mesh->regions()[0];
	}

	size_t materials = _mode & InfoMode::SEPARATE_MATERIALS ? _mesh->materials().size() : 1;

	if (region->elements().size() && region->elements()[0]->params()) {
		std::vector<std::vector<eslocal>*> materialData, bodyData;

		for (size_t r = 0; r < _regions.size(); r++) {
			materialData.push_back(new std::vector<eslocal>());
			bodyData.push_back(new std::vector<eslocal>());
		}

		for (size_t e = 0; e < region->elements().size(); e++) {
			size_t regionOffset = 0, material = -1, body = -1;
			body = region->elements()[e]->param(Element::Params::BODY);
			if (_mode & InfoMode::SEPARATE_BODIES) {
				regionOffset += body * materials;
			}
			material = region->elements()[e]->param(Element::Params::MATERIAL);
			if (_mode & InfoMode::SEPARATE_MATERIALS) {
				regionOffset += material;
			}

			for (auto d = region->elements()[e]->domains().begin(); d != region->elements()[e]->domains().end(); ++d) {
				materialData[regionOffset]->push_back(material);
				bodyData[regionOffset]->push_back(body);
			}
		}

		for (size_t r = 0; r < _regions.size(); r++) {
			_regions[r].data.elementDataInteger["material"] = std::make_pair(1, materialData[r]);
			_regions[r].data.elementDataInteger["body"] = std::make_pair(1, bodyData[r]);
		}
	}

	if (region->settings.size() <= step) {
		return;
	}

	double value;

	for (auto it = region->settings[step].begin(); it != region->settings[step].end(); ++it) {
		const std::vector<Property> &pGroup = _mesh->propertyGroup(it->first);
		if (pGroup.front() != it->first) {
			continue;
		}

		std::vector<std::vector<double>*> rData;

		for (size_t r = 0; r < _regions.size(); r++) {
			rData.push_back(new std::vector<double>());
			rData.back()->reserve(pGroup.size() * _regions[r].elementsTypes.size());
			std::stringstream ss; ss << it->first;
			_regions[r].data.elementDataDouble[ss.str().substr(0, ss.str().find_last_of("_"))] = std::make_pair(pGroup.size(), rData.back());
		}

		for (size_t e = 0; e < region->elements().size(); e++) {
			size_t regionOffset = 0;
			if ((_mode & InfoMode::SEPARATE_BODIES) && region->elements()[e]->params()) {
				regionOffset += region->elements()[e]->param(Element::Params::BODY) * materials;
			}
			if ((_mode & InfoMode::SEPARATE_MATERIALS) && region->elements()[e]->params()) {
				regionOffset += region->elements()[e]->param(Element::Params::MATERIAL);
			}

			for (auto p = pGroup.begin(); p != pGroup.end(); ++p) {
				value = 0;
				for (size_t n = 0; n < region->elements()[e]->nodes(); n++) {
					value += region->elements()[e]->sumProperty(*p, n, step, 0, 0, 0);
				}
				value /= region->elements()[e]->nodes();

				for (auto d = region->elements()[e]->domains().begin(); d != region->elements()[e]->domains().end(); ++d) {
					rData[regionOffset]->push_back(value);
				}
			}
		}
	}
}

void DistributedInfo::addSolution(const std::vector<Solution*> &solution)
{
	for (size_t r = 0; r < _regions.size(); r++) {

		for (size_t s = 0; s < solution.size(); s++) {

			std::vector<double> *rData = new std::vector<double>(solution[s]->properties.size() * _regions[r].coordinates.size() / 3);

			#pragma omp parallel for
			for (size_t d = 0; d < _mesh->parts(); d++) {
				for (size_t i = 0; i < _cIndices[r][d].size(); i++) {
					for (size_t p = 0; p < solution[s]->properties.size(); p++) {
						(*rData)[solution[s]->properties.size() * (i + _cOffset[r][d]) + p] = solution[s]->get(p, d, _mesh->coordinates().localIndex(_cIndices[r][d][i], d));
					}
				}
			}

			switch (solution[s]->eType) {
			case ElementType::ELEMENTS:
				_regions[r].data.elementDataDouble[solution[s]->name] = std::make_pair(solution[s]->properties.size(), rData);
				break;
			case ElementType::NODES:
				_regions[r].data.pointDataDouble[solution[s]->name] = std::make_pair(solution[s]->properties.size(), rData);
				break;
			default:
				ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: cannot store this type solution.";
			}

		}

	}
}


