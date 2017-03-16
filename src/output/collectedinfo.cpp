
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

CollectedInfo::CollectedInfo(const Mesh *mesh)
: RegionInfo(mesh)
{

}

CollectedInfo::CollectedInfo(const Mesh *mesh, size_t body)
: RegionInfo(mesh, body)
{
	prepare(_mesh->elements(), _mesh->getBodies()[body], _mesh->getBodies()[body + 1]);
}

CollectedInfo::CollectedInfo(const Mesh *mesh, const Region* region)
: RegionInfo(mesh, region)
{
	prepare(region->elements(), 0, region->elements().size());
}

RegionInfo* CollectedInfo::deriveRegion(const Region *region) const
{
	return new CollectedInfo(_mesh, region);
}

RegionInfo* CollectedInfo::copyWithoutMesh() const
{
	return new CollectedInfo(_mesh);
}

void CollectedInfo::prepare(const std::vector<Element*> &region, size_t begin, size_t end)
{
	size_t threads = environment->OMP_NUM_THREADS;

	// Collect coordinates
	std::vector<size_t> distribution = Esutils::getDistribution(threads, begin, end);
	std::vector<std::vector<esglobal> > sIDs(threads);
	std::vector<esglobal> rIDs;
	std::vector<double> sCoordinates;
	std::vector<double> rCoordinates;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			for (size_t n = 0; n < region[i]->nodes(); n++) {
				if (_mesh->nodes()[region[i]->node(n)]->clusters().front() == environment->MPIrank) {
					sIDs[t].push_back(_mesh->coordinates().globalIndex(region[i]->node(n)));
				}
			}
		}
	}

	for (size_t t = 1; t < threads; t++) {
		sIDs[0].insert(sIDs[0].end(), sIDs[t].begin(), sIDs[t].end());
	}
	std::sort(sIDs[0].begin(), sIDs[0].end());
	Esutils::removeDuplicity(sIDs[0]);

	rCoordinates.reserve(3 * sIDs[0].size());
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

	std::vector<esglobal> permutation(rIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esglobal i, esglobal j) {
		return rIDs[i] < rIDs[j];
	});

	coordinates.resize(rCoordinates.size());
	distribution = Esutils::getDistribution(threads, rIDs.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			coordinates[3 * i + 0] = rCoordinates[3 * permutation[i] + 0];
			coordinates[3 * i + 1] = rCoordinates[3 * permutation[i] + 1];
			coordinates[3 * i + 2] = rCoordinates[3 * permutation[i] + 2];
		}
	}

	distribution = Esutils::getDistribution(threads, begin, end);
	std::vector<std::vector<esglobal> > sElements(threads);
	std::vector<std::vector<eslocal> > sElementsTypes(threads);
	std::vector<std::vector<eslocal> > sElementsNodes(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			sElementsTypes[t].push_back(region[e]->vtkCode());
			sElementsNodes[t].push_back(region[e]->nodes());
			for (size_t n = 0; n < region[e]->nodes(); n++) {
				sElements[t].push_back(_mesh->coordinates().globalIndex(region[e]->node(n)));
			}
		}
	}

	for (size_t t = 1; t < threads; t++) {
		sElements[0].insert(sElements[0].end(), sElements[t].begin(), sElements[t].end());
		sElementsTypes[0].insert(sElementsTypes[0].end(), sElementsTypes[t].begin(), sElementsTypes[t].end());
		sElementsNodes[0].insert(sElementsNodes[0].end(), sElementsNodes[t].begin(), sElementsNodes[t].end());
	}

	if (!Communication::gatherUnknownSize(sElementsTypes[0], elementsTypes)) {
		ESINFO(ERROR) << "ESPRESO internal error while collecting elements types of a region.";
	}
	if (!Communication::gatherUnknownSize(sElementsNodes[0], elementsNodes)) {
		ESINFO(ERROR) << "ESPRESO internal error while collecting elements nodes of a region.";
	}
	if (!Communication::gatherUnknownSize(sElements[0], elements)) {
		ESINFO(ERROR) << "ESPRESO internal error while collecting elements of a region.";
	}

	distribution = Esutils::getDistribution(threads, elements.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			elements[i] = std::lower_bound(permutation.begin(), permutation.end(), elements[i], [&] (esglobal index, esglobal value) {
				return rIDs[index] < value;
			}) - permutation.begin();
		}
	}

	distribution = Esutils::getDistribution(threads, elementsNodes.size());
	std::vector<eslocal> nodesOffsets(threads);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		eslocal offset = 0;
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			offset += elementsNodes[i];
		}
		nodesOffsets[t] = offset;
	}
	Esutils::sizesToOffsets(nodesOffsets);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		eslocal offset = nodesOffsets[t];
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			elementsNodes[i] += offset;
			offset += elementsNodes[i] - offset;
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
	const Region *region = _region;
	if (region == NULL) {
		region = _mesh->regions()[0];
	}

	if (region->settings.size() <= step) {
		return;
	}
	for (auto it = region->settings[step].begin(); it != region->settings[step].end(); ++it) {
		std::vector<double> sValues;

		sValues.reserve(region->elements().size());
		for (size_t e = 0; e < region->elements().size(); e++) {
			sValues.push_back(0);
			for (size_t n = 0; n < region->elements()[e]->nodes(); n++) {
				sValues.back() += region->elements()[e]->sumProperty(it->first, n, step, 0);
			}
			sValues.back() /= region->elements()[e]->nodes();
			sValues.insert(sValues.end(), region->elements()[e]->domains().size() - 1, sValues.back());
		}

		std::vector<double> *values = new std::vector<double>();
		if (!Communication::gatherUnknownSize(sValues, *values)) {
			ESINFO(ERROR) << "ESPRESO internal error while collecting region data values.";
		}

		std::stringstream ss;
		ss << it->first;
		data.elementDataDouble[ss.str()] = std::make_pair(1, values);
	}
}

void CollectedInfo::addSolution(const std::vector<Solution*> &solution)
{
	if (_body == -1 && _region == NULL) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: cannot store solution.";
	}

	bool status = true;
	for (size_t i = 0; i < solution.size(); i++) {

		const std::vector<Element*> *elements = NULL;
		if (_region != NULL) {
			elements = &_region->elements();
		}
		if (_body != -1) {
			switch (solution[i]->eType) {
			case ElementType::ELEMENTS:
				elements = &_mesh->elements();
				break;
			case ElementType::NODES:
				elements = &_mesh->nodes();
				break;
			default:
				ESINFO(GLOBAL_ERROR) << "Implement storing solution for other element types.";
			}
		}

		size_t threads = environment->OMP_NUM_THREADS;
		std::vector<size_t> distribution = Esutils::getDistribution(threads, elements->size());
		std::vector<double> sData(elements->size() * solution[i]->properties);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

				switch (solution[i]->eType) {
				case ElementType::ELEMENTS:
					for (size_t p = 0; p < solution[i]->properties; p++) {
						sData[e * solution[i]->properties] += solution[i]->get(p, 0, e);
					}
					break;
				case ElementType::NODES:
					for (auto d = (*elements)[e]->domains().begin(); d != (*elements)[e]->domains().end(); ++d) {
						for (size_t p = 0; p < solution[i]->properties; p++) {
							sData[e * solution[i]->properties] += solution[i]->get(p, *d, _mesh->coordinates().localIndex(e, *d));
						}
					}
					break;
				}

			}
		}

		switch (solution[i]->eType) {
		case ElementType::ELEMENTS: {
			std::vector<double> *collected = new std::vector<double>();
			status = status && Communication::gatherUnknownSize(sData, *collected);
			data.elementDataDouble[solution[i]->name] = std::make_pair(solution[i]->properties, collected);
		} break;
		case ElementType::NODES: {
			std::vector<double> collected;
			status = status && Communication::gatherUnknownSize(sData, collected);
			if (!environment->MPIrank) {
				std::vector<double> *averaged = new std::vector<double>();
				data.pointDataDouble[solution[i]->name] = std::make_pair(solution[i]->properties, averaged);
				averaged->resize(coordinates.size() * solution[i]->properties / 3);

				size_t threads = environment->OMP_NUM_THREADS;
				std::vector<size_t> distribution = Esutils::getDistribution(threads, _globalIDsMap.size());
				#pragma omp parallel for
				for (size_t t = 0; t < threads; t++) {
					for (size_t id = distribution[t]; id < distribution[t + 1]; id++) {

						for (size_t p = 0; p < solution[i]->properties; p++) {
							for (esglobal j = (id == 0 ? 0 : _globalIDsMap[id - 1]); j < _globalIDsMap[id]; j++) {
								(*averaged)[id * solution[i]->properties + p] += collected[_globalIDs[j] * solution[i]->properties + p];
							}
							(*averaged)[id * solution[i]->properties + p] /= _globalIDsMultiplicity[id];
						}

					}
				}
			}
		} break;
		default:
			ESINFO(GLOBAL_ERROR) << "Implement storing solution for other element types.";
		}
	}

	if (!status) {
		ESINFO(ERROR) << "ESPRESO internal error while collection solution.";
	}
}


