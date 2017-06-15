
#include "collectedinfo.h"

#include "../configuration/environment.h"

#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"

#include "../basis/point/point.h"
#include "../mesh/elements/element.h"
#include "../mesh/settings/property.h"
#include "../mesh/structures/coordinates.h"
#include "../mesh/structures/mesh.h"
#include "../mesh/structures/region.h"
#include "../mesh/structures/elementtypes.h"

#include "../assembler/solution.h"

#include <numeric>

using namespace espreso::output;

CollectedInfo::CollectedInfo(const Mesh *mesh, InfoMode mode)
: MeshInfo(mesh, mode)
{
	if (_mode | InfoMode::PREPARE) {
		prepare(_mesh->elements());
	}
}

CollectedInfo::CollectedInfo(const Mesh *mesh, const Region* region, InfoMode mode)
: MeshInfo(mesh, region, mode)
{
	if (_mode | InfoMode::PREPARE) {
		prepare(region->elements());
	}
}

MeshInfo* CollectedInfo::deriveRegion(const Region *region) const
{
	return new CollectedInfo(_mesh, region);
}

MeshInfo* CollectedInfo::copyWithoutMesh() const
{
	CollectedInfo* copy = new CollectedInfo(_mesh, _mode & ~MeshInfo::PREPARE);
	copy->_regions.push_back(RegionData());
	return copy;
}

void CollectedInfo::prepare(const std::vector<Element*> &region)
{
	size_t threads = environment->OMP_NUM_THREADS;

	size_t bodies    = _mode & InfoMode::SEPARATE_BODIES    ? _mesh->bodies()           : 1;
	size_t materials = _mode & InfoMode::SEPARATE_MATERIALS ? _mesh->materials().size() : 1;

	_regions.resize(bodies * materials);

	// divide elements to regions and gather them
	std::vector<size_t> distribution = Esutils::getDistribution(threads, region.size());
	std::vector<std::vector<std::vector<esglobal> > > sElements(_regions.size(), std::vector<std::vector<esglobal> >(threads));
	std::vector<std::vector<std::vector<eslocal> > > sElementsTypes(_regions.size(), std::vector<std::vector<esglobal> >(threads));
	std::vector<std::vector<std::vector<eslocal> > > sElementsNodes(_regions.size(), std::vector<std::vector<esglobal> >(threads));

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			size_t regionOffset = 0;
			if ((_mode & InfoMode::SEPARATE_BODIES) && region[e]->params()) {
				regionOffset += region[e]->param(Element::Params::BODY) * materials;
			}
			if ((_mode & InfoMode::SEPARATE_MATERIALS) && region[e]->params()) {
				regionOffset += region[e]->param(Element::Params::MATERIAL);
			}
			sElementsTypes[regionOffset][t].push_back(region[e]->vtkCode());
			sElementsNodes[regionOffset][t].push_back(region[e]->nodes());
			for (size_t n = 0; n < region[e]->nodes(); n++) {
				sElements[regionOffset][t].push_back(_mesh->coordinates().globalIndex(region[e]->node(n)));
			}
		}
	}

	for (size_t r = 0; r < _regions.size(); r++) {
		for (size_t t = 1; t < threads; t++) {
			sElements[r][0].insert(sElements[r][0].end(), sElements[r][t].begin(), sElements[r][t].end());
			sElementsTypes[r][0].insert(sElementsTypes[r][0].end(), sElementsTypes[r][t].begin(), sElementsTypes[r][t].end());
			sElementsNodes[r][0].insert(sElementsNodes[r][0].end(), sElementsNodes[r][t].begin(), sElementsNodes[r][t].end());
		}
	}

	for (size_t r = 0; r < _regions.size(); r++) {
		if (!Communication::gatherUnknownSize(sElementsTypes[r][0], _regions[r].elementsTypes)) {
			ESINFO(ERROR) << "ESPRESO internal error while collecting elements types of a region.";
		}
		if (!Communication::gatherUnknownSize(sElementsNodes[r][0], _regions[r].elementsNodes)) {
			ESINFO(ERROR) << "ESPRESO internal error while collecting elements nodes of a region.";
		}
		if (!Communication::gatherUnknownSize(sElements[r][0], _regions[r].elements)) {
			ESINFO(ERROR) << "ESPRESO internal error while collecting elements of a region.";
		}
	}

	for (size_t r = 0; r < _regions.size(); r++) {
		distribution = Esutils::getDistribution(threads, _regions[r].elementsNodes.size());
		std::vector<eslocal> nodesOffsets(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			eslocal offset = 0;
			for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
				offset += _regions[r].elementsNodes[i];
			}
			nodesOffsets[t] = offset;
		}
		Esutils::sizesToOffsets(nodesOffsets);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			eslocal offset = nodesOffsets[t];
			for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
				_regions[r].elementsNodes[i] += offset;
				offset += _regions[r].elementsNodes[i] - offset;
			}
		}
	}


	// Gather coordinates and divide them to regions
	std::vector<std::vector<esglobal> > sIDs(threads);
	std::vector<esglobal> rIDs;
	std::vector<double> sCoordinates;
	std::vector<double> rCoordinates;

	distribution = Esutils::getDistribution(threads, region.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			for (size_t n = 0; n < region[i]->nodes(); n++) {
				sIDs[t].push_back(_mesh->coordinates().globalIndex(region[i]->node(n)));
			}
		}
	}

	for (size_t t = 1; t < threads; t++) {
		sIDs[0].insert(sIDs[0].end(), sIDs[t].begin(), sIDs[t].end());
	}
	std::sort(sIDs[0].begin(), sIDs[0].end());
	Esutils::removeDuplicity(sIDs[0]);

	sCoordinates.reserve(3 * sIDs[0].size());
	for (size_t i = 0; i < sIDs[0].size(); i++) {
		eslocal n = _mesh->coordinates().clusterIndex(sIDs[0][i]);
		sCoordinates.insert(sCoordinates.end(), { _mesh->coordinates()[n].x, _mesh->coordinates()[n].y, _mesh->coordinates()[n].z });
	}

	if (!Communication::gatherUnknownSize(sIDs[0], rIDs)) {
		ESINFO(ERROR) << "ESPRESO internal error while collecting points global IDs of a region.";
	}
	if (!Communication::gatherUnknownSize(sCoordinates, rCoordinates)) {
		ESINFO(ERROR) << "ESPRESO internal error while collecting points coordinates of a region.";
	}

	std::vector<esglobal> gPermutation(rIDs.size());
	std::iota(gPermutation.begin(), gPermutation.end(), 0);
	std::sort(gPermutation.begin(), gPermutation.end(), [&] (esglobal i, esglobal j) {
		return rIDs[i] < rIDs[j];
	});

	_cIndices.resize(_regions.size());
	for (size_t r = 0; r < _regions.size(); r++) {
		std::vector<std::vector<esglobal> > rNodes(threads);
		distribution = Esutils::getDistribution(threads, _regions[r].elementsTypes.size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {

			size_t begin = distribution[t] ? _regions[r].elementsNodes[distribution[t] - 1] : 0;
			for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

				for (eslocal n = begin; n < _regions[r].elementsNodes[i]; n++) {
					rNodes[t].push_back(_regions[r].elements[n]);
				}
				begin = _regions[r].elementsNodes[i];

			}
		}

		for (size_t t = 1; t < threads; t++) {
			rNodes[0].insert(rNodes[0].end(), rNodes[t].begin(), rNodes[t].end());
		}
		std::sort(rNodes[0].begin(), rNodes[0].end());
		Esutils::removeDuplicity(rNodes[0]);
		_cIndices[r] = rNodes[0];

		_regions[r].coordinates.resize(3 * _cIndices[r].size());
		distribution = Esutils::getDistribution(threads, _cIndices[r].size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
				size_t index = std::lower_bound(gPermutation.begin(), gPermutation.end(), _cIndices[r][i], [&] (esglobal p, esglobal i) { return rIDs[p] < i; }) - gPermutation.begin();
				_regions[r].coordinates[3 * i + 0] = rCoordinates[3 * gPermutation[index] + 0];
				_regions[r].coordinates[3 * i + 1] = rCoordinates[3 * gPermutation[index] + 1];
				_regions[r].coordinates[3 * i + 2] = rCoordinates[3 * gPermutation[index] + 2];
			}
		}

		distribution = Esutils::getDistribution(threads, _regions[r].elements.size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
				_regions[r].elements[i] = std::lower_bound(_cIndices[r].begin(), _cIndices[r].end(), _regions[r].elements[i]) - _cIndices[r].begin();
			}
		}
	}

	// PREPARE SOLUTION VECTOR AVERAGING
	std::vector<int> sMultiplicity(_mesh->coordinates().clusterSize()), rMultiplicity;
	distribution = Esutils::getDistribution(threads, _mesh->coordinates().clusterSize());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			sMultiplicity[i] = _mesh->nodes()[i]->domains().size();
		}
	}

	if (!Communication::gatherUnknownSize(_mesh->coordinates().clusterToGlobal(), rIDs)) {
		ESINFO(ERROR) << "ESPRESO internal error while collecting cluster to global mapping.";
	}
	if (!Communication::gatherUnknownSize(sMultiplicity, rMultiplicity)) {
		ESINFO(ERROR) << "ESPRESO internal error while collecting IDs multiplicity.";
	}

	_globalIDs.resize(rIDs.size());
	std::iota(_globalIDs.begin(), _globalIDs.end(), 0);
	std::sort(_globalIDs.begin(), _globalIDs.end(), [&] (esglobal i, esglobal j) {
		return rIDs[i] < rIDs[j];
	});

	distribution = Esutils::getDistribution(threads, rIDs.size());
	for (size_t t = 1; t < threads; t++) {
		while (distribution[t] < _globalIDs.size() && rIDs[_globalIDs[distribution[t] - 1]] == rIDs[_globalIDs[distribution[t]]]) {
			distribution[t]++;
		}
	}

	std::vector<std::vector<esglobal> > pMap(threads);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		if (distribution[t] < distribution[t + 1]) {
			pMap[t].push_back(distribution[t] + 1);
		}
		for (size_t i = distribution[t] + 1; i < distribution[t + 1]; i++) {
			if (rIDs[_globalIDs[i - 1]] != rIDs[_globalIDs[i]]) {
				pMap[t].insert(pMap[t].end(), rIDs[_globalIDs[i]] - rIDs[_globalIDs[i - 1]] - 1, pMap[t].back());
				pMap[t].push_back(pMap[t].back());
			}
			pMap[t].back()++;
		}
	}

	for (size_t t = 0; t < threads; t++) {
		_globalIDsMap.insert(_globalIDsMap.end(), pMap[t].begin(), pMap[t].end());
	}

	distribution = Esutils::getDistribution(threads, _globalIDsMap.size());
	_globalIDsMultiplicity.resize(_globalIDsMap.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			for (esglobal p = (i == 0 ? 0 : _globalIDsMap[i - 1]); p < _globalIDsMap[i]; p++) {
				_globalIDsMultiplicity[i] += rMultiplicity[_globalIDs[p]];
			}
		}
	}
}

void CollectedInfo::addGeneralInfo()
{

}

void CollectedInfo::addSettings(size_t step)
{
	size_t materials = _mode & InfoMode::SEPARATE_MATERIALS ? _mesh->materials().size() : 1;

	const Region *region = _region;
	if (region == NULL) {
		region = _mesh->regions()[0];
	}

	std::vector<std::vector<eslocal>> sMaterialData(_regions.size()), sBodyData(_regions.size());

	for (size_t e = 0; e < region->elements().size(); e++) {
		size_t regionOffset = 0, material = -1, body = -1;
		if ((_mode & InfoMode::SEPARATE_BODIES) && region->elements()[e]->params()) {
			body = region->elements()[e]->param(Element::Params::BODY);
			regionOffset += body * materials;
		}
		if ((_mode & InfoMode::SEPARATE_MATERIALS) && region->elements()[e]->params()) {
			material = region->elements()[e]->param(Element::Params::MATERIAL);
			regionOffset += material;
		}

		sMaterialData[regionOffset].push_back(material);
		sBodyData[regionOffset].push_back(body);
	}

	for (size_t r = 0; r < _regions.size(); r++) {
		std::vector<eslocal> *values = new std::vector<eslocal>();
		if (!Communication::gatherUnknownSize(sMaterialData[r], *values)) {
			ESINFO(ERROR) << "ESPRESO internal error while collecting region data values.";
		}

		_regions[r].data.elementDataInteger["material"] = std::make_pair(1, values);
	}
	for (size_t r = 0; r < _regions.size(); r++) {
		std::vector<eslocal> *values = new std::vector<eslocal>();
		if (!Communication::gatherUnknownSize(sBodyData[r], *values)) {
			ESINFO(ERROR) << "ESPRESO internal error while collecting region data values.";
		}

		_regions[r].data.elementDataInteger["body"] = std::make_pair(1, values);
	}

	if (region->settings.size() <= step) {
		return;
	}

	for (auto it = region->settings[step].begin(); it != region->settings[step].end(); ++it) {
		const std::vector<Property> &pGroup = _mesh->propertyGroup(it->first);
		if (pGroup.front() != it->first) {
			continue;
		}

		std::vector<std::vector<double> > sValues(_regions.size());

		for (size_t e = 0; e < region->elements().size(); e++) {
			size_t regionOffset = 0;
			if ((_mode & InfoMode::SEPARATE_BODIES) && region->elements()[e]->params()) {
				regionOffset += region->elements()[e]->param(Element::Params::BODY) * materials;
			}
			if ((_mode & InfoMode::SEPARATE_MATERIALS) && region->elements()[e]->params()) {
				regionOffset += region->elements()[e]->param(Element::Params::MATERIAL);
			}

			for (auto p = pGroup.begin(); p != pGroup.end(); ++p) {
				sValues[regionOffset].push_back(0);
				for (size_t n = 0; n < region->elements()[e]->nodes(); n++) {
					sValues[regionOffset].back() += region->elements()[e]->sumProperty(*p, n, step, 0, 0, 0);
				}
				sValues[regionOffset].back() /= region->elements()[e]->nodes();
			}
		}

		for (size_t r = 0; r < _regions.size(); r++) {
			std::vector<double> *values = new std::vector<double>();
			if (!Communication::gatherUnknownSize(sValues[r], *values)) {
				ESINFO(ERROR) << "ESPRESO internal error while collecting region data values.";
			}

			std::stringstream ss;
			ss << it->first;
			_regions[r].data.elementDataDouble[ss.str().substr(0, ss.str().find_last_of("_"))] = std::make_pair(pGroup.size(), values);
		}
	}
}

void CollectedInfo::addSolution(const std::vector<Solution*> &solution)
{
	if (_region != NULL) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: cannot store solution.";
	}

	bool status = true;
	size_t threads = environment->OMP_NUM_THREADS;

	size_t materials = _mode & InfoMode::SEPARATE_MATERIALS ? _mesh->materials().size() : 1;

	for (size_t i = 0; i < solution.size(); i++) {


		if (solution[i]->eType == ElementType::ELEMENTS) {
			std::vector<std::vector<double> > sData(_regions.size());

			std::vector<std::vector<std::vector<double> > > rData(_regions.size(), std::vector<std::vector<double> >(threads));
			std::vector<size_t> distribution = Esutils::getDistribution(threads, _mesh->elements().size());

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
					size_t regionOffset = 0;
					if ((_mode & InfoMode::SEPARATE_BODIES) && _mesh->elements()[e]->params()) {
						regionOffset += _mesh->elements()[e]->param(Element::Params::BODY) * materials;
					}
					if ((_mode & InfoMode::SEPARATE_MATERIALS) && _mesh->elements()[e]->params()) {
						regionOffset += _mesh->elements()[e]->param(Element::Params::MATERIAL);
					}

					eslocal d = std::lower_bound(_mesh->getPartition().begin(), _mesh->getPartition().end(), e + 1) - _mesh->getPartition().begin() - 1;
					for (size_t p = 0; p < solution[i]->properties.size(); p++) {
						rData[regionOffset][t].push_back(solution[i]->get(p, d, e - _mesh->getPartition()[d]));
					}

				}
			}

			for (size_t r = 0; r < _regions.size(); r++) {
				for (size_t t = 0; t < threads; t++) {
					sData[r].insert(sData[r].end(), rData[r][t].begin(), rData[r][t].end());
				}

				std::vector<double> *collected = new std::vector<double>();
				status = status && Communication::gatherUnknownSize(sData[r], *collected);
				_regions[r].data.elementDataDouble[solution[i]->name] = std::make_pair(solution[i]->properties.size(), collected);
			}
		}

		if (solution[i]->eType == ElementType::NODES) {
			std::vector<double> sData(solution[i]->properties.size() * _mesh->nodes().size());

			std::vector<size_t> distribution = Esutils::getDistribution(threads, _mesh->nodes().size());
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {

					for (auto d = _mesh->nodes()[n]->domains().begin(); d != _mesh->nodes()[n]->domains().end(); ++d) {
						for (size_t p = 0; p < solution[i]->properties.size(); p++) {
							sData[n * solution[i]->properties.size() + p] += solution[i]->get(p, *d, _mesh->coordinates().localIndex(n, *d));
						}
					}
				}
			}

			std::vector<double> collected;
			status = status && Communication::gatherUnknownSize(sData, collected);

			for (size_t r = 0; r < _regions.size(); r++) {
				std::vector<double> *averaged = new std::vector<double>();
				_regions[r].data.pointDataDouble[solution[i]->name] = std::make_pair(solution[i]->properties.size(), averaged);
				averaged->resize(_cIndices[r].size() * solution[i]->properties.size());

				std::vector<size_t> distribution = Esutils::getDistribution(threads, _cIndices[r].size());
				#pragma omp parallel for
				for (size_t t = 0; t < threads; t++) {
					for (size_t id = distribution[t]; id < distribution[t + 1]; id++) {

						for (size_t p = 0; p < solution[i]->properties.size(); p++) {
							for (esglobal j = _cIndices[r][id] ? _globalIDsMap[_cIndices[r][id] - 1] : 0; j < _globalIDsMap[_cIndices[r][id]]; j++) {
								(*averaged)[id * solution[i]->properties.size() + p] += collected[_globalIDs[j] * solution[i]->properties.size() + p];
							}
							(*averaged)[id * solution[i]->properties.size() + p] /= _globalIDsMultiplicity[_cIndices[r][id]];
						}

					}
				}
			}
		}
	}

	if (!status) {
		ESINFO(ERROR) << "ESPRESO internal error while collection solution.";
	}
}


