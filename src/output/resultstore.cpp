
#include "resultstore.h"

#include "collectedinfo.h"
#include "distributedinfo.h"

#include <numeric>

#include "../configuration/output.h"

#include "../basis/point/point.h"
#include "../mesh/elements/element.h"
#include "../mesh/elements/point/node.h"
#include "../mesh/elements/line/line2.h"
#include "../mesh/structures/coordinates.h"
#include "../mesh/structures/mesh.h"
#include "../mesh/structures/region.h"
#include "../mesh/settings/property.h"

#include "../assembler/step.h"
#include "../assembler/solution.h"
#include "../assembler/instance.h"

#include "../solver/generic/SparseMatrix.h"

#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"

using namespace espreso::output;

ResultStore::ResultStore(const OutputConfiguration &output, const Mesh *mesh, const std::string &path, MeshInfo::InfoMode mode)
: Store(output), _mesh(mesh), _path(path), _meshInfo(NULL)
{
	if (_configuration.separate_bodies) {
		mode = mode | MeshInfo::SEPARATE_BODIES;
	}
	if (_configuration.separate_materials) {
		mode = mode | MeshInfo::SEPARATE_MATERIALS;
	}
	if (mode & MeshInfo::InfoMode::PREPARE) {
		if (_configuration.collected) {
			_meshInfo = new CollectedInfo(_mesh, mode);
		} else {
			_meshInfo = new DistributedInfo(_mesh, _configuration.domain_shrink_ratio, _configuration.cluster_shrink_ratio, mode);
		}
	}
}

ResultStore::~ResultStore()
{
	if (_meshInfo != NULL) {
		delete _meshInfo;
	}
}

std::vector<std::string> ResultStore::store(const std::string &name, const Step &step, const MeshInfo *meshInfo)
{
	std::string root = Esutils::createDirectory({ "results", "step" + std::to_string(step.step), "substep" + std::to_string(step.substep) });
	std::vector<std::string> files;

	if (meshInfo->distributed()) {
		std::string prefix = Esutils::createDirectory({ "results", "step" + std::to_string(step.step), "substep" + std::to_string(step.substep), std::to_string(environment->MPIrank) });
		for (size_t r = 0; r < meshInfo->regions(); r++) {
			store(prefix + name + std::to_string(r), meshInfo->region(r));
			if (!environment->MPIrank) {
				files.push_back(linkClusters(root, name + std::to_string(r), meshInfo->region(r)));
			}
		}

	} else {
		if (!environment->MPIrank) {
			for (size_t r = 0; r < meshInfo->regions(); r++) {
				files.push_back(store(root + name + std::to_string(r), meshInfo->region(r)));
			}
		}
	}

	return files;
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

void ResultStore::storeSettings(const std::vector<size_t> &steps)
{
	Step step;
	std::vector<std::string> files;

	for (size_t i = 0; i < steps.size(); i++) {
		step.step = steps[i];

		_meshInfo->addSettings(i);
		files = store("mesh", step, _meshInfo);
		_meshInfo->clearData();
		_settings.push_back(std::make_pair(step, files));
	}

	MeshInfo *region;
	for (size_t r = 2; r < _mesh->regions().size(); r++) {
		region = _meshInfo->deriveRegion(_mesh->regions()[r]);
		for (size_t i = 0; i < steps.size(); i++) {
			step.step = steps[i];

			region->addSettings(i);
			files = store(_mesh->regions()[r]->name, step, region);
			region->clearData();
			_settings.push_back(std::make_pair(step, files));
		}
		delete region;
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
		linkSteps("solution", _solutions);
	}
	_solutions.clear();
}

void ResultStore::storeSolution(const Step &step, const std::vector<Solution*> &solution)
{
	_meshInfo->addSolution(solution);
	std::vector<std::string> files = store("solution", step, _meshInfo);
	_meshInfo->clearData();

	_solutions.push_back(std::make_pair(step, files));
}

void ResultStore::finalize()
{
	if (!environment->MPIrank) {
		if (_solutions.size()) {
			linkSteps("solution", _solutions);
		}
		if (_settings.size()) {
			linkSteps("settings", _settings);
		}
		if (_FETIdata.size()) {
			linkSteps("feti_data", _FETIdata);
		}
	}
}


void ResultStore::storeFETIData(const Step &step, const Instance &instance)
{
	storeElementInfo(step);
	storeFixPoints(step);
	storeCorners(step);
	return; // FIX STORING B1
	storeDirichlet(step, instance);
	storeLambdas(step, instance);
}

void ResultStore::storeElementInfo(const Step &step)
{
	_meshInfo->addGeneralInfo();
	std::vector<std::string> files = store("mesh_info", step, _meshInfo);
	_meshInfo->clearData();
	_FETIdata.push_back(std::make_pair(step, files));
}

void ResultStore::storeFixPoints(const Step &step)
{
	std::vector<Element*> fixPoints;
	for (size_t p = 0; p < _mesh->fixPoints().size(); p++) {
		fixPoints.insert(fixPoints.end(), _mesh->fixPoints(p).begin(), _mesh->fixPoints(p).end());
	}

	std::sort(fixPoints.begin(), fixPoints.end());
	Esutils::removeDuplicity(fixPoints);

	Region region(fixPoints);

	MeshInfo *info = _meshInfo->deriveRegion(&region);
	std::vector<std::string> files = store("fix_points", step, info);
	delete info;
	_FETIdata.push_back(std::make_pair(step, files));
}

void ResultStore::storeCorners(const Step &step)
{
	std::vector<Element*> corners = _mesh->corners();

	Region region(corners);

	MeshInfo *info = _meshInfo->deriveRegion(&region);
	std::vector<std::string> files = store("corners", step, info);
	delete info;
	_FETIdata.push_back(std::make_pair(step, files));
}

void ResultStore::storeDirichlet(const Step &step, const Instance &instance)
{
	MeshInfo *info = _meshInfo->copyWithoutMesh();
	for (size_t p = 0; p < instance.properties.size(); p++) {
		std::vector<double> *values = new std::vector<double>();
		Point point;

		for (size_t d = 0; d < instance.domains; d++) {
			for (size_t i = 0; i < instance.B1[d].I_row_indices.size() && instance.B1[d].I_row_indices[i] <= instance.block[Instance::CONSTRAINT::DIRICHLET]; i++) {
				const Element *e = _mesh->getDOFsElement(d, instance.B1[d].J_col_indices[i] - 1);
				point = _meshInfo->shrink(_mesh->coordinates()[e->node(0)], d);
				info->_regions[0].coordinates.insert(info->_regions[0].coordinates.end(), { point.x, point.y, point.z });
				values->push_back(instance.B1c[d][i]);
			}
		}

		info->_regions[0].elementsTypes.resize(values->size(), NodeVTKCode);
		info->_regions[0].elementsNodes.resize(values->size());
		info->_regions[0].elements.resize(values->size());
		std::iota(info->_regions[0].elementsNodes.begin(), info->_regions[0].elementsNodes.end(), 1);
		std::iota(info->_regions[0].elements.begin(), info->_regions[0].elements.end(), 0);
		std::stringstream ss;
		ss << instance.properties[p];
		info->_regions[0].data.pointDataDouble[ss.str()] = std::make_pair(1, values);
		store("DIRICHLET_" + ss.str(), step, info);
		info->clearData();
	}
	delete info;
}

void ResultStore::storeLambdas(const Step &step, const Instance &instance)
{
	if (!_meshInfo->distributed()) {
		// it is pointless to store lambdas for collected result
		return;
	}
	MeshInfo *info = _meshInfo->copyWithoutMesh();

	std::vector<int> neighbours(environment->MPIsize);
	std::iota(neighbours.begin(), neighbours.end(), 0);

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

	std::vector<double> duplicity(lMap.size()), *values = new std::vector<double>(lMap.size());

	info->_regions[0].coordinates.resize(6 * lMap.size());
	info->_regions[0].elementsTypes.resize(lMap.size(), -1);
	info->_regions[0].elementsNodes.reserve(lMap.size());
	info->_regions[0].elements.reserve(2 * lMap.size());

	for (size_t i = 0; i < lMap.size(); i++) {
		info->_regions[0].elementsNodes.push_back(2 * i + 2);
		info->_regions[0].elements.push_back(2 * i);
		info->_regions[0].elements.push_back(2 * i + 1);
	}

	auto lIndex = [&] (esglobal lambda) {
		std::pair<esglobal, esglobal> pair(lambda, 0);
		auto it = std::lower_bound(lMap.begin(), lMap.end(), pair);
		if (it == lMap.end() || it->first != lambda) {
			ESINFO(ERROR) << "Invalid lamdas in B1.";
		}
		return it - lMap.begin();
	};

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
			point = info->shrink(_mesh->coordinates()[e->node(0)], d);
			size_t index = lIndex(it->front());
			for (size_t c = 1; c < it->size(); c++) {
				if ((*it)[c] == environment->MPIrank) {
					if (info->_regions[0].elementsTypes[index] != -1) {
						index++;
						if (lMap[index].first != it->front()) {
							ESINFO(ERROR) << "Broken B1.";
						}
					}
					info->_regions[0].elementsTypes[index] = Line2VTKCode;
					info->_regions[0].coordinates[6 * index + 0] = point.x;
					info->_regions[0].coordinates[6 * index + 1] = point.y;
					info->_regions[0].coordinates[6 * index + 2] = point.z;
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
				info->_regions[0].coordinates[6 * lindex + 3] = rPoints[(*it)[c]][index].x;
				info->_regions[0].coordinates[6 * lindex + 4] = rPoints[(*it)[c]][index].y;
				info->_regions[0].coordinates[6 * lindex + 5] = rPoints[(*it)[c]][index].z;
			}
		}
	}

	for (size_t i = 0; i < lMap.size(); i++) {
		if (i + 1 < lMap.size() && lMap[i].first == lMap[i + 1].first) {
			Point a(info->_regions[0].coordinates[6 * (i + 0) + 0], info->_regions[0].coordinates[6 * (i + 0) + 1], info->_regions[0].coordinates[6 * (i + 0) + 2]);
			Point b(info->_regions[0].coordinates[6 * (i + 1) + 0], info->_regions[0].coordinates[6 * (i + 1) + 1], info->_regions[0].coordinates[6 * (i + 1) + 2]);
			Point length = b - a;
			Point ax = a + length * duplicity[i];
			Point bx = b - length * duplicity[i + 1];
			info->_regions[0].coordinates[6 * (i + 0) + 3] = ax.x;
			info->_regions[0].coordinates[6 * (i + 0) + 4] = ax.y;
			info->_regions[0].coordinates[6 * (i + 0) + 5] = ax.z;
			info->_regions[0].coordinates[6 * (i + 1) + 3] = bx.x;
			info->_regions[0].coordinates[6 * (i + 1) + 4] = bx.y;
			info->_regions[0].coordinates[6 * (i + 1) + 5] = bx.z;
			i++;
		} else {
			Point a(info->_regions[0].coordinates[6 * (i + 0) + 0], info->_regions[0].coordinates[6 * (i + 0) + 1], info->_regions[0].coordinates[6 * (i + 0) + 2]);
			Point b(info->_regions[0].coordinates[6 * (i + 0) + 3], info->_regions[0].coordinates[6 * (i + 0) + 4], info->_regions[0].coordinates[6 * (i + 0) + 5]);
			Point length = b - a;
			Point ax = a + length * duplicity[i];
			info->_regions[0].coordinates[6 * (i + 0) + 3] = ax.x;
			info->_regions[0].coordinates[6 * (i + 0) + 4] = ax.y;
			info->_regions[0].coordinates[6 * (i + 0) + 5] = ax.z;
		}
	}

	info->_regions[0].data.elementDataDouble["values"] = std::make_pair(1, values);
	store("EQUALITY_CONSTRAINTS", step, info);
	delete info;
}






