
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

static size_t getDataSize(const espreso::Mesh &mesh, const std::vector<std::vector<double> > &data, espreso::ElementType eType)
{
	if (data.size() != mesh.parts()) {
		ESINFO(ERROR) << "ESPRESO internal error: invalid statistic data - number of domains not match.";
	}

	size_t dataSize;
	for (size_t p = 0; p < mesh.parts(); p++) {
		switch (eType) {
			case ElementType::NODES:
				if (data[p].size() % mesh.coordinates().localSize(p) != 0) {
					ESINFO(ERROR) << "ESPRESO internal error: invalid statistic data - size of domain data is not multiplication of local size.";
				}
				if (p && dataSize != data[p].size() / mesh.coordinates().localSize(p)) {
					ESINFO(ERROR) << "ESPRESO internal error: invalid statistic data - size of domain data is not the same for all domains.";
				}
				dataSize = data[p].size() / mesh.coordinates().localSize(p);
				break;
			case ElementType::ELEMENTS:
				if (data[p].size() % (mesh.getPartition()[p + 1] - mesh.getPartition()[p]) != 0) {
					ESINFO(ERROR) << "ESPRESO internal error: invalid statistic data - size of domain data is not multiplication of local size.";
				}
				if (p && dataSize != data[p].size() / (mesh.getPartition()[p + 1] - mesh.getPartition()[p])) {
					ESINFO(ERROR) << "ESPRESO internal error: invalid statistic data - size of domain data is not the same for all domains.";
				}
				dataSize = data[p].size() / (mesh.getPartition()[p + 1] - mesh.getPartition()[p]);
				break;
			default:
				break;
		}

	}

	return dataSize;
}

Statistic::Statistic(StatisticalData statistics, Operation operation, ElementType eType, const Mesh &mesh, const std::vector<std::vector<double> > &data, const std::vector<size_t> &offsets, const std::vector<Region*> &selection)
: _statistics(statistics), _operation(operation), _eType(eType), _computed(false), _mesh(mesh), _data(data), _offsets(offsets), _selection(selection)
{
	_dataSize = getDataSize(mesh, data, eType);
}

Statistic::Statistic(ElementType eType, const Mesh &mesh, const std::vector<std::vector<double> > &data, size_t dataSize)
: _statistics(StatisticalData::MIN | StatisticalData::MAX | StatisticalData::AVERAGE | StatisticalData::NORM), _operation(Operation::AVERAGE), _eType(eType), _computed(false), _mesh(mesh), _data(data)
{
	for (size_t r = 0; r < _mesh.monitoredRegions().size(); r++) {
		_selection.push_back(_mesh.monitoredRegions()[r]);
	}

	_dataSize = dataSize;
	_offsets.resize(_dataSize);
	std::iota(_offsets.begin(), _offsets.end(), 0);
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

	if (!(statistics & _statistics)) {
		ESINFO(ERROR) << "ESPRESO internal error: request for not pre-computed statistics.";
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
	_computed = true;
	if (!_selection.size()) {
		return;
	}

	if (_eType != ElementType::NODES) {
		ESINFO(GLOBAL_ERROR) << "Implement monitoring non-node elements.";
	}

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh.neighbours().begin(), _mesh.neighbours().end(), neighbour) - _mesh.neighbours().begin();
	};

	std::vector<bool> _regions(_mesh.regions().size(), false);
	for (size_t s = 0; s < _selection.size(); s++) {
		// mark selected regions
		size_t index = std::find(_mesh.regions().begin(), _mesh.regions().end(), _selection[s]) - _mesh.regions().begin();
		_regions[index] = true;
	}

	std::vector<Element*> _elements;
	const std::vector<Element*> &elements = _regions[1] ? _mesh.nodes() : _regions[0] ? _mesh.elements() : _elements;

	auto c2l = [ & ] (eslocal domain, size_t e) -> eslocal {
		switch (_eType) {
			case ElementType::ELEMENTS:
				return e - _mesh.getPartition()[domain];
			case ElementType::NODES:
				return _mesh.coordinates().localIndex(elements[e]->node(0), domain);
				break;
			default:
				ESINFO(GLOBAL_ERROR) << "Implement statistic for new element types.";
				break;
		}
		return -1;
	};

	if (!_regions[0] && !_regions[1]) {
		for (size_t r = 0; r < _mesh.regions().size(); r++) {
			if (_regions[r]) {
				if (_mesh.regions()[r]->eType != ElementType::NODES) {
					for (size_t e = 0; e < _mesh.regions()[r]->elements().size(); e++) {
						for (size_t n = 0; n < _mesh.regions()[r]->elements()[e]->nodes(); n++) {
							_elements.push_back(_mesh.nodes()[_mesh.regions()[r]->elements()[e]->node(n)]);
						}
					}
				} else {
					_elements.insert(_elements.end(), _mesh.regions()[r]->elements().begin(), _mesh.regions()[r]->elements().end());
				}
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

			if (_eType == ElementType::NODES && elements[e]->parentElements().size() == 0) {
				// mesh generator can generate dangling nodes -> skip them
				continue;
			}
			if (elements[e]->clusters().size() == 1) {
				continue;
			}

			if (elements[e]->clusters().front() == environment->MPIrank) {
				// process only element that have to be send to lower rank
				for (size_t c = 1; c < elements[e]->clusters().size(); c++) {
					rIDs[t][n2i(elements[e]->clusters()[c])].push_back(e);
					rOffset[n2i(elements[e]->clusters()[c])][t]++;
				}
				continue;
			}

			size_t n = n2i(elements[e]->clusters().front());


			for (size_t i = 0; i < _offsets.size(); i++) {
				sBuffer[t][n].push_back(0);
				for (auto d = elements[e]->domains().begin(); d != elements[e]->domains().end(); ++d) {
					sBuffer[t][n].back() += _data[*d][_dataSize * c2l(*d, e) + _offsets[i]];
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
			for (size_t j = 0; j < _offsets.size(); j++) {
				permuted.push_back(rBuffer[n][_offsets.size() * permutation[i] + j]);
			}
		}
		rBuffer[n].swap(permuted);
	}

	// MIN, MAX, AVERAGE, NORM, SQUARES, ELEMENTS
	std::vector<double> initData = { std::numeric_limits<double>::max(), std::numeric_limits<double>::min(), 0, 0, 0, 0 };

	// thread x region x offsets x data
	std::vector<std::vector<std::vector<std::vector<double> > > > tData(threads, std::vector<std::vector<std::vector<double> > >(_selection.size(), std::vector<std::vector<double> >(_offsets.size(), initData)));

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<double> value(_offsets.size());

		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			if (_eType == ElementType::NODES && elements[e]->parentElements().size() == 0) {
				// mesh generator can generate dangling nodes -> skip them
				continue;
			}

			if (elements[e]->clusters().front() != environment->MPIrank) {
				continue;
			}

			std::fill(value.begin(), value.end(), 0);
			for (auto d = elements[e]->domains().begin(); d != elements[e]->domains().end(); ++d) {
				for (size_t i = 0; i < _offsets.size(); i++) {
					value[i] += _data[*d][_dataSize * c2l(*d, e) + _offsets[i]];
				}
			}

			for (size_t c = 1; c < elements[e]->clusters().size(); c++) {
				size_t n = n2i(elements[e]->clusters()[c]);
				for (size_t i = 0; i < _offsets.size(); i++) {
					value[i] += rBuffer[n][_dataSize * rOffset[n][t] + i];
				}
				rOffset[n][t]++;
			}

			if (_operation == Operation::AVERAGE) {
				for (size_t i = 0; i < _offsets.size(); i++) {
					value[i] /= elements[e]->numberOfGlobalDomains();
				}
			}

			for (size_t i = 0; i < _offsets.size(); i++) {
				for (size_t r = 0; r < _selection.size(); r++) {
					if (std::binary_search(elements[e]->regions().begin(), elements[e]->regions().end(), _selection[r])) {
						if (_statistics & StatisticalData::MIN) {
							tData[t][r][i][0] = std::min(tData[t][r][i][0], value[i]);
						}

						if (_statistics & StatisticalData::MAX) {
							tData[t][r][i][1] = std::max(tData[t][r][i][1], value[i]);
						}

						if (_statistics & StatisticalData::AVERAGE) {
							tData[t][r][i][2] += value[i];
						}

						if (_statistics & (StatisticalData::NORM || StatisticalData::SQUARES)) {
							tData[t][r][i][3] += value[i] * value[i];
						}

						tData[t][r][i][5] += 1;
					}
				}
			}

		}
	}

	std::vector<double> cData;
	cData.reserve(_selection.size() * _offsets.size() * initData.size());
	_results.resize(_selection.size());
	for (size_t r = 0; r < _selection.size(); r++) {
		_results[r].resize(_offsets.size());
		for (size_t i = 0; i < _offsets.size(); i++) {
			cData.insert(cData.end(), initData.begin(), initData.end());
			_results[r][i] = initData;
		}
	}

	for (size_t t = 0; t < threads; t++) {
		for (size_t i = 0, offset = 0; i < _offsets.size(); i++) {
			for (size_t r = 0; r < _selection.size(); r++, offset++) {
				cData[initData.size() * offset + 0] = std::min(cData[initData.size() * offset + 0], tData[t][r][i][0]);
				cData[initData.size() * offset + 1] = std::max(cData[initData.size() * offset + 1], tData[t][r][i][1]);
				cData[initData.size() * offset + 2] += tData[t][r][i][2];
				cData[initData.size() * offset + 3] += tData[t][r][i][3];
				cData[initData.size() * offset + 4] += tData[t][r][i][3];
				cData[initData.size() * offset + 5] += tData[t][r][i][5];
			}
		}
	}

	std::vector<double> gData;
	if (!Communication::gatherUnknownSize(cData, gData)) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error while gathering global statistic.";
	}

	for (int r = 0; r < environment->MPIsize && !environment->MPIrank; r++) {
		for (size_t i = 0, offset = 0; i < _offsets.size(); i++) {
			for (size_t s = 0; s < _selection.size(); s++, offset++) {
				_results[s][i][0] = std::min(_results[s][i][0], gData[r * cData.size() + initData.size() * offset + 0]);
				_results[s][i][1] = std::max(_results[s][i][1], gData[r * cData.size() + initData.size() * offset + 1]);
				_results[s][i][2] += gData[r * cData.size() + initData.size() * offset + 2];
				_results[s][i][3] += gData[r * cData.size() + initData.size() * offset + 3];
				_results[s][i][4] += gData[r * cData.size() + initData.size() * offset + 3];
				_results[s][i][5] += gData[r * cData.size() + initData.size() * offset + 5];
			}
		}
	}

	if (!environment->MPIrank) {
		for (size_t i = 0, offset = 0; i < _offsets.size(); i++) {
			for (size_t s = 0; s < _selection.size(); s++, offset++) {
				_results[s][i][2] /= _results[s][i][5];
				_results[s][i][3] = std::sqrt(_results[s][i][3]);
			}
		}
	}

	for (size_t r = 0; r < _selection.size(); r++) {
		for (size_t i = 0; i < _offsets.size(); i++) {
			MPI_Bcast(_results[r][i].data(), _results[r][i].size() * sizeof(double), MPI_BYTE, 0, environment->MPICommunicator);
		}
	}
}




