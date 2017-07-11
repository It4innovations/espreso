
#include "statistic.h"

#include "../basis/logging/logging.h"
#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"

#include "../configuration/environment.h"

#include "../mesh/structures/coordinates.h"
#include "../mesh/structures/mesh.h"
#include "../mesh/structures/region.h"
#include "../mesh/structures/elementtypes.h"
#include "../mesh/elements/element.h"

#include <limits>
#include <numeric>

using namespace espreso;


Statistic::Statistic(ElementType eType, const Mesh &mesh, const std::vector<std::vector<double> > &data, const std::vector<Property> &properties)
: _operation(Operation::AVERAGE),
  _eType(eType), _dataSize(properties.size()), _computed(false), _mesh(mesh), _data(data)
{
	for (size_t r = 0; r < _mesh.monitoredRegions().size(); r++) {
		if (eType == ElementType::ELEMENTS) {
			if (_mesh.monitoredRegions()[r]->eType == ElementType::ELEMENTS) {
				_selection.push_back(_mesh.monitoredRegions()[r]);
			}
		} else {
			_selection.push_back(_mesh.monitoredRegions()[r]);
		}
	}
}

double Statistic::get(const Region* region, size_t offset, StatisticalData statistics)
{
	if (!_computed) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: request for not pre-computed statistic. Call compute statistical data first!";
	}
	size_t index = std::find(_selection.begin(), _selection.end(), region) - _selection.begin();
	if (index == _selection.size() || !StringCompare::caseInsensitiveEq(_selection[index]->name, region->name)) {
		ESINFO(ERROR) << "Cannot compute statistic for region '" << region->name << "'. Elements of the region do not have the requested property.";
	}

	switch (statistics) {
	case StatisticalData::MIN    : return _results[index][offset][0];
	case StatisticalData::MAX    : return _results[index][offset][1];
	case StatisticalData::AVERAGE: return _results[index][offset][2];
	case StatisticalData::NORM   : return _results[index][offset][3];
	case StatisticalData::SQUARES: return _results[index][offset][4];
	default:
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid statistics";
		return 0;
	}
}

void Statistic::compute()
{
	if (_computed) {
		return;
	}
	_computed = true;
	if (!_selection.size()) {
		return;
	}

	switch (_eType) {
	case ElementType::NODES:
		computeNodes();
		break;
	case ElementType::ELEMENTS:
		computeElements();
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Implement statistic for new element types.";
	}
}

void Statistic::computeNodes()
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh.neighbours().begin(), _mesh.neighbours().end(), neighbour) - _mesh.neighbours().begin();
	};

	std::vector<Element*> _elements;
	const std::vector<Element*> &elements = std::find(_selection.begin(), _selection.end(), _mesh.regions()[1]) != _selection.end() ? _mesh.nodes() : _elements;

	if (std::find(_selection.begin(), _selection.end(), _mesh.regions()[1]) == _selection.end()) {
		for (size_t s = 0; s < _selection.size(); s++) {
			if (_selection[s]->eType != ElementType::NODES) {
				for (size_t e = 0; e < _selection[s]->elements().size(); e++) {
					for (size_t n = 0; n < _selection[s]->elements()[e]->nodes(); n++) {
						_elements.push_back(_mesh.nodes()[_selection[s]->elements()[e]->node(n)]);
					}
				}
			} else {
				_elements.insert(_elements.end(), _selection[s]->elements().begin(), _selection[s]->elements().end());
			}
		}
		std::sort(_elements.begin(), _elements.end());
		Esutils::removeDuplicity(_elements);
	}

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());

	// thread x neighbour x data
	std::vector<std::vector<std::vector<double> > > sBuffer(threads, std::vector<std::vector<double> >(_mesh.neighbours().size()));
	std::vector<std::vector<std::vector<esglobal> > > rIDs(threads, std::vector<std::vector<esglobal> >(_mesh.neighbours().size()));
	std::vector<std::vector<size_t> > rOffset(_mesh.neighbours().size(), std::vector<size_t>(threads));

	std::vector<std::vector<double> > rBuffer(_mesh.neighbours().size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			if (elements[e]->parentElements().size() == 0) {
				// mesh generator can generate dangling nodes -> skip them
				continue;
			}
			if (elements[e]->clusters().size() == 1) {
				// will be counted later
				continue;
			}

			if (elements[e]->clusters().front() == environment->MPIrank) {
				// process only element that have to be send to lower rank
				for (size_t c = 1; c < elements[e]->clusters().size(); c++) {
					rIDs[t][n2i(elements[e]->clusters()[c])].push_back(e);
					rOffset[n2i(elements[e]->clusters()[c])][t]++;
				}
			} else {
				size_t n = n2i(elements[e]->clusters().front());
				for (size_t i = 0; i < _dataSize; i++) {
					sBuffer[t][n].push_back(0);
					for (auto d = elements[e]->domains().begin(); d != elements[e]->domains().end(); ++d) {
						sBuffer[t][n].back() += _data[*d][_dataSize * _mesh.coordinates().localIndex(elements[e]->node(0), *d) + i];
					}
				}
			}

		}
	}

	#pragma omp parallel for
	for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
		for (size_t t = 1; t < threads; t++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
			rIDs[0][n].insert(rIDs[0][n].end(), rIDs[t][n].begin(), rIDs[t][n].end());
		}
		Esutils::sizesToOffsets(rOffset[n]);
	}

	if (!Communication::receiveUpperUnknownSize(sBuffer[0], rBuffer, _mesh.neighbours())) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error in computing statistic.";
	}

	#pragma omp parallel for
	for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
		std::vector<eslocal> permutation(rIDs[0][n].size());
		std::iota(permutation.begin(), permutation.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
			return elements[rIDs[0][n][i]]->clusterOffset(_mesh.neighbours()[n]) < elements[rIDs[0][n][j]]->clusterOffset(_mesh.neighbours()[n]);
		});

		std::vector<double> permuted;
		permuted.reserve(permutation.size());
		for (size_t i = 0; i < permutation.size(); i++) {
			for (size_t j = 0; j < _dataSize; j++) {
				permuted.push_back(rBuffer[n][_dataSize * permutation[i] + j]);
			}
		}
		rBuffer[n].swap(permuted);
	}

	// MIN, MAX, AVERAGE, NORM, SQUARES, ELEMENTS
	std::vector<double> initData = { std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), 0, 0, 0, 0 };

	// thread x region x offsets x data
	std::vector<std::vector<std::vector<std::vector<double> > > > tData(threads, std::vector<std::vector<std::vector<double> > >(_selection.size(), std::vector<std::vector<double> >(_dataSize + 1, initData)));

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<double> value(_dataSize + 1);
		std::vector<Region*> eregions;

		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			if (elements[e]->parentElements().size() == 0) {
				// mesh generator can generate dangling nodes -> skip them
				continue;
			}

			if (elements[e]->clusters().front() != environment->MPIrank) {
				continue;
			}

			std::fill(value.begin(), value.end(), 0);
			for (auto d = elements[e]->domains().begin(); d != elements[e]->domains().end(); ++d) {
				for (size_t i = 0; i < _dataSize; i++) {
					value[i] += _data[*d][_dataSize * _mesh.coordinates().localIndex(elements[e]->node(0), *d) + i];
				}
			}

			for (size_t c = 1; c < elements[e]->clusters().size(); c++) {
				size_t n = n2i(elements[e]->clusters()[c]);
				for (size_t i = 0; i < _dataSize; i++) {
					value[i] += rBuffer[n][_dataSize * rOffset[n][t] + i];
				}
				rOffset[n][t]++;
			}

			if (_operation == Operation::AVERAGE) {
				for (size_t i = 0; i < _dataSize; i++) {
					value[i] /= elements[e]->numberOfGlobalDomains();
				}
			}

			for (size_t i = 0; i < _dataSize; i++) {
				value[_dataSize] += value[i] * value[i];
			}
			value[_dataSize] = std::sqrt(value[_dataSize]);

			eregions.clear();
			eregions.insert(eregions.end(), elements[e]->regions().begin(), elements[e]->regions().end());
			std::for_each(elements[e]->parentElements().begin(), elements[e]->parentElements().end(), [&] (Element *parent) {
				eregions.insert(eregions.end(), parent->regions().begin(), parent->regions().end());
			});
			std::for_each(elements[e]->parentFaces().begin(), elements[e]->parentFaces().end(), [&] (Element *parent) {
				eregions.insert(eregions.end(), parent->regions().begin(), parent->regions().end());
			});
			std::for_each(elements[e]->parentEdges().begin(), elements[e]->parentEdges().end(), [&] (Element *parent) {
				eregions.insert(eregions.end(), parent->regions().begin(), parent->regions().end());
			});
			std::sort(eregions.begin(), eregions.end());
			Esutils::removeDuplicity(eregions);
			for (size_t i = 0; i <= _dataSize; i++) {
				for (size_t r = 0; r < _selection.size(); r++) {
					if (std::binary_search(eregions.begin(), eregions.end(), _selection[r])) {
						tData[t][r][i][0] = std::min(tData[t][r][i][0], value[i]);
						tData[t][r][i][1] = std::max(tData[t][r][i][1], value[i]);
						tData[t][r][i][2] += value[i];
						tData[t][r][i][3] += value[i] * value[i];
						tData[t][r][i][4] += value[i] * value[i];
						tData[t][r][i][5] += 1;
					}
				}
			}

		}
	}

	std::vector<double> cData;
	cData.reserve(_selection.size() * (_dataSize + 1) * initData.size());
	_results.resize(_selection.size());
	for (size_t r = 0; r < _selection.size(); r++) {
		_results[r].resize(_dataSize + 1);
		for (size_t i = 0; i <= _dataSize; i++) {
			cData.insert(cData.end(), initData.begin(), initData.end());
			_results[r][i] = initData;
		}
	}

	for (size_t t = 0; t < threads; t++) {
		for (size_t i = 0, offset = 0; i <= _dataSize; i++) {
			for (size_t r = 0; r < _selection.size(); r++, offset++) {
				cData[initData.size() * offset + 0] = std::min(cData[initData.size() * offset + 0], tData[t][r][i][0]);
				cData[initData.size() * offset + 1] = std::max(cData[initData.size() * offset + 1], tData[t][r][i][1]);
				cData[initData.size() * offset + 2] += tData[t][r][i][2];
				cData[initData.size() * offset + 3] += tData[t][r][i][3];
				cData[initData.size() * offset + 4] += tData[t][r][i][4];
				cData[initData.size() * offset + 5] += tData[t][r][i][5];
			}
		}
	}

	std::vector<double> gData;
	if (!Communication::gatherUnknownSize(cData, gData)) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error while gathering global statistic.";
	}

	for (int r = 0; r < environment->MPIsize && !environment->MPIrank; r++) {
		for (size_t i = 0, offset = 0; i <= _dataSize; i++) {
			for (size_t s = 0; s < _selection.size(); s++, offset++) {
				_results[s][i][0] = std::min(_results[s][i][0], gData[r * cData.size() + initData.size() * offset + 0]);
				_results[s][i][1] = std::max(_results[s][i][1], gData[r * cData.size() + initData.size() * offset + 1]);
				_results[s][i][2] += gData[r * cData.size() + initData.size() * offset + 2];
				_results[s][i][3] += gData[r * cData.size() + initData.size() * offset + 3];
				_results[s][i][4] += gData[r * cData.size() + initData.size() * offset + 4];
				_results[s][i][5] += gData[r * cData.size() + initData.size() * offset + 5];
			}
		}
	}

	if (!environment->MPIrank) {
		for (size_t i = 0; i <= _dataSize; i++) {
			for (size_t s = 0; s < _selection.size(); s++) {
				_results[s][i][2] /= _results[s][i][5];
				_results[s][i][3] = std::sqrt(_results[s][i][3]);
			}
		}
	}

	for (size_t r = 0; r < _selection.size(); r++) {
		for (size_t i = 0; i <= _dataSize; i++) {
			MPI_Bcast(_results[r][i].data(), _results[r][i].size() * sizeof(double), MPI_BYTE, 0, environment->MPICommunicator);
		}
	}
}

void Statistic::computeElements()
{
	// MIN, MAX, AVERAGE, NORM, SQUARES, ELEMENTS
	std::vector<double> initData = { std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), 0, 0, 0, 0 };

	// part x region x offsets x data
	std::vector<std::vector<std::vector<std::vector<double> > > > tData(_mesh.parts(), std::vector<std::vector<std::vector<double> > >(_selection.size(), std::vector<std::vector<double> >(_dataSize + 1, initData)));

	#pragma omp parallel for
	for (size_t p = 0; p < _mesh.parts(); p++) {
		std::vector<double> value(_dataSize + 1);

		for (size_t e = _mesh.getPartition()[p]; e < _mesh.getPartition()[p + 1]; e++) {

			std::fill(value.begin(), value.end(), 0);
			for (size_t i = 0; i < _dataSize; i++) {
				value[i] = _data[p][_dataSize * (e - _mesh.getPartition()[p]) + i];
			}

			for (size_t i = 0; i < _dataSize; i++) {
				value[_dataSize] += value[i] * value[i];
			}
			value[_dataSize] = std::sqrt(value[_dataSize]);

			for (size_t i = 0; i <= _dataSize; i++) {
				for (size_t r = 0; r < _selection.size(); r++) {
					if (std::binary_search(_mesh.elements()[e]->regions().begin(), _mesh.elements()[e]->regions().end(), _selection[r])) {
						tData[p][r][i][0] = std::min(tData[p][r][i][0], value[i]);
						tData[p][r][i][1] = std::max(tData[p][r][i][1], value[i]);
						tData[p][r][i][2] += value[i];
						tData[p][r][i][3] += value[i] * value[i];
						tData[p][r][i][4] += value[i] * value[i];

						tData[p][r][i][5] += 1;
					}
				}
			}

		}
	}

	std::vector<double> cData;
	cData.reserve(_selection.size() * (_dataSize + 1) * initData.size());

	_results.resize(_selection.size());
	for (size_t r = 0; r < _selection.size(); r++) {
		_results[r].resize(_dataSize + 1);
		for (size_t i = 0; i <= _dataSize; i++) {
			cData.insert(cData.end(), initData.begin(), initData.end());
			_results[r][i] = initData;
		}
	}

	for (size_t p = 0; p < _mesh.parts(); p++) {
		for (size_t i = 0, offset = 0; i <= _dataSize; i++) {
			for (size_t r = 0; r < _selection.size(); r++, offset++) {
				cData[initData.size() * offset + 0] = std::min(cData[initData.size() * offset + 0], tData[p][r][i][0]);
				cData[initData.size() * offset + 1] = std::max(cData[initData.size() * offset + 1], tData[p][r][i][1]);
				cData[initData.size() * offset + 2] += tData[p][r][i][2];
				cData[initData.size() * offset + 3] += tData[p][r][i][3];
				cData[initData.size() * offset + 4] += tData[p][r][i][4];
				cData[initData.size() * offset + 5] += tData[p][r][i][5];
			}
		}
	}

	std::vector<double> gData;
	if (!Communication::gatherUnknownSize(cData, gData)) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error while gathering global statistic.";
	}


	for (int r = 0; r < environment->MPIsize && !environment->MPIrank; r++) {
		for (size_t i = 0, offset = 0; i <= _dataSize; i++) {
			for (size_t s = 0; s < _selection.size(); s++, offset++) {
				_results[s][i][0] = std::min(_results[s][i][0], gData[r * cData.size() + initData.size() * offset + 0]);
				_results[s][i][1] = std::max(_results[s][i][1], gData[r * cData.size() + initData.size() * offset + 1]);
				_results[s][i][2] += gData[r * cData.size() + initData.size() * offset + 2];
				_results[s][i][3] += gData[r * cData.size() + initData.size() * offset + 3];
				_results[s][i][4] += gData[r * cData.size() + initData.size() * offset + 4];
				_results[s][i][5] += gData[r * cData.size() + initData.size() * offset + 5];
			}
		}
	}

	if (!environment->MPIrank) {
		for (size_t i = 0; i <= _dataSize; i++) {
			for (size_t s = 0; s < _selection.size(); s++) {
				_results[s][i][2] /= _results[s][i][5];
				_results[s][i][3] = std::sqrt(_results[s][i][3]);
			}
		}
	}

	for (size_t r = 0; r < _selection.size(); r++) {
		for (size_t i = 0; i <= _dataSize; i++) {
			MPI_Bcast(_results[r][i].data(), _results[r][i].size() * sizeof(double), MPI_BYTE, 0, environment->MPICommunicator);
		}
	}
}




