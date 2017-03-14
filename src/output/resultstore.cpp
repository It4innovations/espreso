
#include "resultstore.h"

#include <numeric>

#include "../configuration/output.h"
#include "../configuration/environment.h"

#include "../basis/point/point.h"
#include "../mesh/elements/element.h"
#include "../mesh/elements/point/node.h"
#include "../mesh/elements/line/line2.h"
#include "../mesh/structures/coordinates.h"
#include "../mesh/structures/mesh.h"
#include "../mesh/structures/region.h"
#include "../mesh/settings/property.h"
#include "../mesh/structures/elementtypes.h"

#include "../assembler/step.h"
#include "../assembler/solution.h"
#include "../assembler/instance.h"

#include "../solver/generic/SparseMatrix.h"

#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"

using namespace espreso::output;

DataArrays::~DataArrays()
{
	for (auto it = elementDataDouble.begin(); it != elementDataDouble.end(); ++it) {
		delete it->second.second;
	}
	for (auto it = elementDataInteger.begin(); it != elementDataInteger.end(); ++it) {
		delete it->second.second;
	}
	for (auto it = pointDataDouble.begin(); it != pointDataDouble.end(); ++it) {
		delete it->second.second;
	}
	for (auto it = pointDataInteger.begin(); it != pointDataInteger.end(); ++it) {
		delete it->second.second;
	}
}

ResultStore::ResultStore(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: Store(output), _mesh(mesh), _path(path)
{

}

ResultStore::~ResultStore()
{

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

void ResultStore::coordinatePreprocessing(const std::vector<std::vector<eslocal> > &indices, std::vector<double> &coordinates, std::vector<size_t> &offsets)
{
	// compute coordinates
	std::vector<std::vector<double> > dCoordinates(_mesh->parts());
	#pragma omp parallel for
	for (size_t p = 0; p < indices.size(); p++) {
		Point point;
		dCoordinates[p].reserve(indices[p].size());
		for (size_t i = 0; i < indices[p].size(); i++) {
			point = _mesh->coordinates()[indices[p][i]];
			point = _clusterCenter + (point - _clusterCenter) * _configuration.cluster_shrink_ratio;
			point = _domainsCenters[p] + (point - _domainsCenters[p]) * _configuration.domain_shrink_ratio;
			dCoordinates[p].insert(dCoordinates[p].end(), { point.x, point.y, point.z });
		}
	}

	coordinates.clear();
	offsets = { 0 };
	for (size_t p = 0; p < _mesh->parts(); p++) {
		coordinates.insert(coordinates.end(), dCoordinates[p].begin(), dCoordinates[p].end());
		offsets.push_back(coordinates.size() / 3);
	}
}

void ResultStore::elementsPreprocessing(const std::vector<Element*> &region, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements)
{
	if (!region.size()) {
		// TODO: collected output
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
	coordinatePreprocessing(rCoordinates, coordinates, offsets);

	for (size_t e = 0, offset = 0; e < region.size(); e++) {
		elementsTypes.insert(elementsTypes.end(), region[e]->domains().size(), region[e]->vtkCode());
		for (auto d = region[e]->domains().begin(); d != region[e]->domains().end(); ++d, offset += region[e]->nodes()) {
			elementsNodes.push_back(offset + region[e]->nodes());
			for (size_t n = 0; n < region[e]->nodes(); n++) {
				eslocal oIndex = std::lower_bound(rCoordinates[*d].begin(), rCoordinates[*d].end(), region[e]->node(n)) - rCoordinates[*d].begin();
				elements.push_back(oIndex + offsets[*d]);
			}
		}
	}
}

void ResultStore::regionData(size_t step, const espreso::Region &region, DataArrays &data)
{
	if (region.settings.size() <= step) {
		return;
	}
	for (auto it = region.settings[step].begin(); it != region.settings[step].end(); ++it) {
		std::vector<double> *values = new std::vector<double>();
		values->reserve(region.elements().size());
		for (size_t e = 0; e < region.elements().size(); e++) {
			values->push_back(0);
			for (size_t n = 0; n < region.elements()[e]->nodes(); n++) {
				values->back() += region.elements()[e]->sumProperty(it->first, n, step, 0);
			}
			values->back() /= region.elements()[e]->nodes();
			values->insert(values->end(), region.elements()[e]->domains().size() - 1, values->back());
		}
		std::stringstream ss;
		ss << it->first;
		data.elementDataDouble[ss.str()] = std::make_pair(1, values);
	}
}

void ResultStore::preprocessing()
{
	if (_configuration.collected) {
		// Implement collected printing
	} else {
		_clusterCenter = computeClusterCenter(_mesh);
		_domainsCenters = computeDomainsCenters(_mesh);
		std::vector<size_t> domainOffset;

		coordinatePreprocessing(_mesh->coordinates().localToCluster(), _coordinates, domainOffset);

		// compute elements
		std::vector<std::vector<eslocal> > dElementsTypes(_mesh->parts());
		std::vector<std::vector<eslocal> > dElementsNodes(_mesh->parts());
		std::vector<std::vector<eslocal> > dElements(_mesh->parts());
		std::vector<size_t> offsets(_mesh->parts());
		#pragma omp parallel for
		for (size_t p = 0; p < _mesh->parts(); p++) {
			size_t offset = 0;
			for (size_t e = _mesh->getPartition()[p]; e < _mesh->getPartition()[p + 1]; e++) {
				dElementsTypes[p].push_back(_mesh->elements()[e]->vtkCode());
				dElementsNodes[p].push_back(_mesh->elements()[e]->nodes());
				offset += _mesh->elements()[e]->nodes();
				for (size_t n = 0; n < _mesh->elements()[e]->nodes(); n++) {
					dElements[p].push_back(_mesh->coordinates().localIndex(_mesh->elements()[e]->node(n), p) + domainOffset[p]);
				}
			}
			offsets[p] = offset;
		}

		Esutils::sizesToOffsets(offsets);
		#pragma omp parallel for
		for (size_t p = 0; p < _mesh->parts(); p++) {
			for (size_t e = _mesh->getPartition()[p], i = 0; e < _mesh->getPartition()[p + 1]; e++, i++) {
				dElementsNodes[p][i] = offsets[p] += dElementsNodes[p][i];
			}
		}

		for (size_t p = 0; p < _mesh->parts(); p++) {
			_elementsTypes.insert(_elementsTypes.end(), dElementsTypes[p].begin(), dElementsTypes[p].end());
			_elementsNodes.insert(_elementsNodes.end(), dElementsNodes[p].begin(), dElementsNodes[p].end());
			_elements.insert(_elements.end(), dElements[p].begin(), dElements[p].end());
		}
	}
}

static void fillMeshSettings(DataArrays &data, const std::vector<double> &coordinates, const espreso::Mesh *mesh)
{
	std::vector<eslocal> *pointIDcluster = new std::vector<eslocal>(coordinates.size() / 3);
	std::vector<eslocal> *pointIDglobal = new std::vector<eslocal>();
	std::vector<eslocal> *elementID = new std::vector<eslocal>(mesh->getPartition().back());
	std::vector<eslocal> *decomposition = new std::vector<eslocal>();

	std::iota(pointIDcluster->begin(), pointIDcluster->end(), 0);
	std::iota(elementID->begin(), elementID->end(), 0);
	for (size_t p = 0; p < mesh->getPartition().size() - 1; p++) {
		decomposition->insert(decomposition->end(), mesh->getPartition()[p + 1] - mesh->getPartition()[p], p);
	}

	pointIDglobal->reserve(coordinates.size() / 3);
	for (size_t p = 0; p < mesh->parts(); p++) {
		for (size_t i = 0; i < mesh->coordinates().localSize(p); i++) {
			pointIDglobal->push_back(mesh->coordinates().globalIndex(i, p));
		}
	}

	data.pointDataInteger["pointIDcluster"] = std::make_pair(1, pointIDcluster);
	data.pointDataInteger["pointIDglobal"] = std::make_pair(1, pointIDglobal);
	data.elementDataInteger["elementID"] = std::make_pair(1, elementID);
	data.elementDataInteger["decomposition"] = std::make_pair(1, decomposition);

	for (int r = 0; r < espreso::environment->MPIsize; r++) {
		std::vector<eslocal> *cluster = new std::vector<eslocal>(coordinates.size() / 3);
		data.pointDataInteger["cluster" + std::to_string(r)] = std::make_pair(1, cluster);
	}

	std::vector<size_t> offsets(mesh->parts());
	for (size_t p = 1; p < mesh->parts(); p++) {
		offsets[p] = offsets[p - 1] + mesh->coordinates().localSize(p - 1);
	}

	for (size_t n = 0; n < mesh->nodes().size(); n++) {
		for (auto c = mesh->nodes()[n]->clusters().begin(); c != mesh->nodes()[n]->clusters().end(); ++c) {
			for (auto d = mesh->nodes()[n]->domains().begin(); d != mesh->nodes()[n]->domains().end(); ++d) {
				(*data.pointDataInteger["cluster" + std::to_string(*c)].second)[offsets[*d] + mesh->coordinates().localIndex(n, *d)] = 1;
			}
		}
	}
}

void ResultStore::storeSettings(const Step &step)
{
	storeSettings(std::vector<size_t>{ step.step });
}

void ResultStore::storeSettings(size_t steps)
{
	std::vector<size_t> _steps(steps);
	std::iota(_steps.begin(), _steps.end(), 0);
	storeSettings(_steps);
}

template <class TData>
std::string ResultStore::store(const std::string &name, const Step &step, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, const TData &data)
{
	std::string root, prefix;
	if (!environment->MPIrank) {
		root = Esutils::createDirectory({ "results", "step" + std::to_string(step.step), "substep" + std::to_string(step.substep) });
	}
	prefix = Esutils::createDirectory({ "results", "step" + std::to_string(step.step), "substep" + std::to_string(step.substep), std::to_string(environment->MPIrank) });

	store(prefix + name, coordinates, elementsTypes, elementsNodes, elements, data);
	if (!environment->MPIrank) {
		linkClusters(root, name, data);
	}
	return root;
}

void ResultStore::storeSettings(const std::vector<size_t> &steps)
{
	Step step;

	DataArrays data;
	fillMeshSettings(data, _coordinates, _mesh);
	for (size_t i = 0; i < steps.size(); i++) {
		step.step = steps[i];
		store("mesh", step, _coordinates, _elementsTypes, _elementsNodes, _elements, data);
	}

	for (size_t r = 0; r < _mesh->regions().size(); r++) {
		std::vector<double> coordinates;
		std::vector<eslocal> elementsTypes, elementsNodes, elements;
		elementsPreprocessing(_mesh->regions()[r]->elements(), coordinates, elementsTypes, elementsNodes, elements);

		for (size_t i = 0; i < steps.size(); i++) {
			DataArrays rData;
			regionData(steps[step.step], *_mesh->regions()[r], rData);

			step.step = steps[i];
			store(_mesh->regions()[r]->name, step, coordinates, elementsTypes, elementsNodes, elements, rData);
		}
	}
}

void ResultStore::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType)
{
	Step step;
	std::vector<Solution*> solution;
	std::vector<Property> props;
	if (dimension == 1) {
		props.push_back(Property::TEMPERATURE);
	}
	if (dimension == 3) {
		props.push_back(Property::DISPLACEMENT_X);
		props.push_back(Property::DISPLACEMENT_Y);
		props.push_back(Property::DISPLACEMENT_Z);
	}

	solution.push_back(new Solution(name, eType, props, values));
	storeSolution(step, solution);
	delete solution.back();
	if (!environment->MPIrank) {
		linkSteps("solution", _steps);
	}
	_steps.clear();
}

void ResultStore::storeSolution(const Step &step, const std::vector<Solution*> &solution)
{
	std::string file;

	file = store("solution", step, _coordinates, _elementsTypes, _elementsNodes, _elements, solution);
	_steps.push_back(std::make_pair(file, step));
}

void ResultStore::finalize()
{
	if (!environment->MPIrank && _steps.size()) {
		linkSteps("solution", _steps);
	}
}


void ResultStore::storeFETIData(const Step &step, const Instance &instance)
{
	storeFixPoints(step);
	storeCorners(step);
	storeDirichlet(step, instance);
	storeLambdas(step, instance);
}

void ResultStore::storeFixPoints(const Step &step)
{
	std::vector<Element*> fixPoints;
	for (size_t p = 0; p < _mesh->fixPoints().size(); p++) {
		fixPoints.insert(fixPoints.end(), _mesh->fixPoints(p).begin(), _mesh->fixPoints(p).end());
	}

	std::sort(fixPoints.begin(), fixPoints.end());
	Esutils::removeDuplicity(fixPoints);

	DataArrays data;
	std::vector<double> coordinates;
	std::vector<eslocal> elementsTypes, elementsNodes, elements;
	elementsPreprocessing(fixPoints, coordinates, elementsTypes, elementsNodes, elements);
	store("fix_points", step, coordinates, elementsTypes, elementsNodes, elements, data);
}

void ResultStore::storeCorners(const Step &step)
{
	std::vector<Element*> corners = _mesh->corners();

	DataArrays data;
	std::vector<double> coordinates;
	std::vector<eslocal> elementsTypes, elementsNodes, elements;
	elementsPreprocessing(corners, coordinates, elementsTypes, elementsNodes, elements);
	store("corners", step, coordinates, elementsTypes, elementsNodes, elements, data);
}

void ResultStore::storeDirichlet(const Step &step, const Instance &instance)
{
	for (size_t p = 0; p < instance.properties.size(); p++) {
		DataArrays data;
		std::vector<double> coordinates, *values = new std::vector<double>();
		std::vector<eslocal> elementsTypes, elementsNodes, elements;
		Point point;

		for (size_t d = 0; d < instance.domains; d++) {
			for (size_t i = 0; instance.B1[d].I_row_indices[i] <= instance.block[Instance::CONSTRAINT::DIRICHLET]; i++) {
				const Element *e = _mesh->getDOFsElement(d, instance.B1[d].J_col_indices[i] - 1);
				point = _mesh->coordinates()[e->node(0)];
				point = _clusterCenter + (point - _clusterCenter) * _configuration.cluster_shrink_ratio;
				point = _domainsCenters[d] + (point - _domainsCenters[d]) * _configuration.domain_shrink_ratio;
				coordinates.insert(coordinates.end(), { point.x, point.y, point.z });
				values->push_back(instance.B1c[d][i]);
			}
		}

		elementsTypes.resize(values->size(), NodeVTKCode);
		elementsNodes.resize(values->size());
		elements.resize(values->size());
		std::iota(elementsNodes.begin(), elementsNodes.end(), 1);
		std::iota(elements.begin(), elements.end(), 0);
		std::stringstream ss;
		ss << instance.properties[p];
		data.pointDataDouble[ss.str()] = std::make_pair(1, values);
		store("DIRICHLET_" + ss.str(), step, coordinates, elementsTypes, elementsNodes, elements, data);
	}
}

void ResultStore::storeLambdas(const Step &step, const Instance &instance)
{
	std::vector<int> neighbours(environment->MPIsize);
	std::iota(neighbours.begin(), neighbours.end(), 0);

	DataArrays data;
	std::vector<std::pair<esglobal, esglobal> > lMap;

	for (size_t i = 0; i < instance.B1clustersMap.size(); i++) {
		if (instance.B1clustersMap[i].front() < instance.block[Instance::CONSTRAINT::DIRICHLET]) {
			continue;
		}
		if (instance.B1clustersMap[i].size() == 2) {
			lMap.push_back(std::make_pair((esglobal)instance.B1clustersMap[i][0], (esglobal)(lMap.size())));
			lMap.push_back(std::make_pair((esglobal)instance.B1clustersMap[i][0], (esglobal)(lMap.size())));
		} else {
			lMap.push_back(std::make_pair((esglobal)instance.B1clustersMap[i][0], (esglobal)(lMap.size())));
		}
	}

	std::vector<double> coordinates(6 * lMap.size()), duplicity(lMap.size()), *values = new std::vector<double>(lMap.size());
	std::vector<eslocal> elementsTypes(lMap.size(), -1), elementsNodes, elements;
	elementsNodes.reserve(lMap.size());
	elements.reserve(2 * lMap.size());

	for (size_t i = 0; i < lMap.size(); i++) {
		elementsNodes.push_back(2 * i + 2);
		elements.push_back(2 * i);
		elements.push_back(2 * i + 1);
	}

	auto lIndex = [&] (esglobal lambda) {
		std::pair<esglobal, esglobal> pair(lambda, 0);
		auto it = std::lower_bound(lMap.begin(), lMap.end(), pair);
		if (it == lMap.end() || it->first != lambda) {
			ESINFO(ERROR) << "Invalid lamdas in B1.";
		}
		return it - lMap.begin();
	};

	for (size_t p = 0; p < instance.properties.size(); p++) {
		std::vector<std::vector<esglobal> > sLambdas(environment->MPIsize);
		std::vector<std::vector<Point> > sPoints(environment->MPIsize);
		std::vector<std::vector<esglobal> > rLambdas(environment->MPIsize);
		std::vector<std::vector<Point> > rPoints(environment->MPIsize);
		Point point;

		for (size_t d = 0; d < instance.domains; d++) {
			auto start = std::upper_bound(instance.B1[d].I_row_indices.begin(), instance.B1[d].I_row_indices.end(), instance.block[Instance::CONSTRAINT::DIRICHLET]);
			auto end = std::upper_bound(instance.B1[d].I_row_indices.begin(), instance.B1[d].I_row_indices.end(), instance.block[Instance::CONSTRAINT::DIRICHLET] + instance.block[Instance::CONSTRAINT::EQUALITY_CONSTRAINTS]);
			for (size_t i = start - instance.B1[d].I_row_indices.begin(); i < end - instance.B1[d].I_row_indices.begin(); i++) {
				auto it = std::lower_bound(instance.B1clustersMap.begin(), instance.B1clustersMap.end(), instance.B1[d].I_row_indices[i] - 1, [&] (const std::vector<esglobal> &v, esglobal i) {
					return v[0] < i;
				});
				const Element *e = _mesh->getDOFsElement(d, instance.B1[d].J_col_indices[i] - 1);
				point = _mesh->coordinates()[e->node(0)];
				point = _clusterCenter + (point - _clusterCenter) * _configuration.cluster_shrink_ratio;
				point = _domainsCenters[d] + (point - _domainsCenters[d]) * _configuration.domain_shrink_ratio;
				size_t index = lIndex(it->front());
				for (size_t c = 1; c < it->size(); c++) {
					if ((*it)[c] == environment->MPIrank) {
						if (elementsTypes[index] != -1) {
							index++;
							if (lMap[index].first != it->front()) {
								ESINFO(ERROR) << "Broken B1.";
							}
						}
						elementsTypes[index] = Line2VTKCode;
						coordinates[6 * index + 0] = point.x;
						coordinates[6 * index + 1] = point.y;
						coordinates[6 * index + 2] = point.z;
						(*values)[index] = instance.B1[d].V_values[i];
						duplicity[index] = instance.B1duplicity[d][i];
					} else {
						sLambdas[(*it)[c]].push_back(it->front());
						sPoints[(*it)[c]].push_back(point);
					}
				}
			}
		}


		if (!Communication::exchangeUnknownSize(sLambdas, rLambdas, neighbours)) {
			ESINFO(ERROR) << "problem while exchanging Lambdas in storeLambdas.";
		}
		if (!Communication::exchangeUnknownSize(sPoints, rPoints, neighbours)) {
			ESINFO(ERROR) << "problem while exchanging Points in storeLambdas.";
		}
		for (int i = 0; i < environment->MPIsize; i++) {
			if (sLambdas[i].size() != rLambdas[i].size() || sPoints[i].size() != rPoints[i].size()) {
				ESINFO(ERROR) << "Lambda indices do not match.";
			}
		}

		for (size_t d = 0; d < instance.domains; d++) {
			auto start = std::upper_bound(instance.B1[d].I_row_indices.begin(), instance.B1[d].I_row_indices.end(), instance.block[Instance::CONSTRAINT::DIRICHLET]);
			auto end = std::upper_bound(instance.B1[d].I_row_indices.begin(), instance.B1[d].I_row_indices.end(), instance.block[Instance::CONSTRAINT::DIRICHLET] + instance.block[Instance::CONSTRAINT::EQUALITY_CONSTRAINTS]);
			for (size_t i = start - instance.B1[d].I_row_indices.begin(); i < end - instance.B1[d].I_row_indices.begin(); i++) {
				auto it = std::lower_bound(instance.B1clustersMap.begin(), instance.B1clustersMap.end(), instance.B1[d].I_row_indices[i] - 1, [&] (const std::vector<esglobal> &v, esglobal i) {
					return v[0] < i;
				});
				size_t lindex = lIndex(it->front());
				for (size_t c = 2; c < it->size(); c++) {
					size_t index = std::find(rLambdas[(*it)[c]].begin(), rLambdas[(*it)[c]].end(), it->front()) - rLambdas[(*it)[c]].begin();
					if (index == rLambdas[(*it)[c]].size() || rLambdas[(*it)[c]][index] != it->front()) {
						ESINFO(ERROR) << "Different Lambdas on neighbour clusters.";
					}
					coordinates[6 * lindex + 3] = rPoints[(*it)[c]][index].x;
					coordinates[6 * lindex + 4] = rPoints[(*it)[c]][index].y;
					coordinates[6 * lindex + 5] = rPoints[(*it)[c]][index].z;
				}
			}
		}

		for (size_t i = 0; i < lMap.size(); i++) {
			if (i + 1 < lMap.size() && lMap[i].first == lMap[i + 1].first) {
				Point a(coordinates[6 * (i + 0) + 0], coordinates[6 * (i + 0) + 1], coordinates[6 * (i + 0) + 2]);
				Point b(coordinates[6 * (i + 1) + 0], coordinates[6 * (i + 1) + 1], coordinates[6 * (i + 1) + 2]);
				Point length = b - a;
				Point ax = a + length * duplicity[i];
				Point bx = b - length * duplicity[i + 1];
				coordinates[6 * (i + 0) + 3] = ax.x;
				coordinates[6 * (i + 0) + 4] = ax.y;
				coordinates[6 * (i + 0) + 5] = ax.z;
				coordinates[6 * (i + 1) + 3] = bx.x;
				coordinates[6 * (i + 1) + 4] = bx.y;
				coordinates[6 * (i + 1) + 5] = bx.z;
				i++;
			} else {
				Point a(coordinates[6 * (i + 0) + 0], coordinates[6 * (i + 0) + 1], coordinates[6 * (i + 0) + 2]);
				Point b(coordinates[6 * (i + 0) + 3], coordinates[6 * (i + 0) + 4], coordinates[6 * (i + 0) + 5]);
				Point length = b - a;
				Point ax = a + length * duplicity[i];
				coordinates[6 * (i + 0) + 3] = ax.x;
				coordinates[6 * (i + 0) + 4] = ax.y;
				coordinates[6 * (i + 0) + 5] = ax.z;
			}
		}

		std::stringstream ss;
		ss << instance.properties[p];
		data.elementDataDouble["values"] = std::make_pair(1, values);
		store("EQUALITY_CONSTRAINTS_" + ss.str(), step, coordinates, elementsTypes, elementsNodes, elements, data);
	}
}






