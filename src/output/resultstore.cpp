
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
#include "../mesh/structures/elementtypes.h"
#include "../mesh/structures/region.h"
#include "../mesh/settings/property.h"

#include "../assembler/step.h"
#include "../assembler/solution.h"
#include "../assembler/instance.h"

#include "../solver/generic/SparseMatrix.h"

#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"

using namespace espreso;

ResultStore::ResultStore(const OutputConfiguration &output, const Mesh *mesh, MeshInfo::InfoMode mode)
: Store(output), _mesh(mesh), _meshInfo(NULL)
{
	if (_configuration.separate_bodies) {
		mode = mode | MeshInfo::SEPARATE_BODIES;
	}
	if (_configuration.separate_materials) {
		mode = mode | MeshInfo::SEPARATE_MATERIALS;
	}
	if (mode & MeshInfo::InfoMode::PREPARE) {
		prepare();
	}
	_mode = mode | MeshInfo::InfoMode::PREPARE;
}

ResultStore::~ResultStore()
{
	if (_meshInfo != NULL) {
		delete _meshInfo;
	}
}

void ResultStore::prepare()
{
	if (_meshInfo == NULL) {
		if (_configuration.collected) {
			_meshInfo = new CollectedInfo(_mesh, _mode);
		} else {
			_meshInfo = new DistributedInfo(_mesh, _configuration.domain_shrink_ratio, _configuration.cluster_shrink_ratio, _mode);
		}
	}
}

std::vector<std::string> ResultStore::store(const std::string &name, const Step &step, const MeshInfo *meshInfo)
{
	std::string root;
	if (_configuration.subsolution) {
		root = Esutils::createDirectory({ Logging::outputRoot(), "PRE_POST_DATA", "step" + std::to_string(step.step), "substep" + std::to_string(step.substep), "iteration" + std::to_string(step.iteration)});
	} else {
		root = Esutils::createDirectory({ Logging::outputRoot(), "PRE_POST_DATA", "step" + std::to_string(step.step), "substep" + std::to_string(step.substep) });
	}
	std::vector<std::string> files;

	if (meshInfo->distributed()) {
		std::string prefix;
		if (_configuration.subsolution) {
			prefix = Esutils::createDirectory({ Logging::outputRoot(), "PRE_POST_DATA", "step" + std::to_string(step.step), "substep" + std::to_string(step.substep), "iteration" + std::to_string(step.iteration), std::to_string(environment->MPIrank) });
		} else {
			prefix = Esutils::createDirectory({ Logging::outputRoot(), "PRE_POST_DATA", "step" + std::to_string(step.step), "substep" + std::to_string(step.substep), std::to_string(environment->MPIrank) });
		}
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

void ResultStore::storeSettings(size_t steps)
{
	if (!_configuration.settings) {
		return;
	}

	prepare();
	Step step;
	std::vector<std::string> files;

	for (size_t i = 0; i < steps; i++) {
		step.step = i;
		step.currentTime = i;

		_meshInfo->addSettings(i);
		files = store("mesh", step, _meshInfo);
		_meshInfo->clearData();
		_settings.push_back(std::make_pair(step, files));
	}

	MeshInfo *region;
	for (size_t r = 2; r < _mesh->regions().size(); r++) {
		region = _meshInfo->deriveRegion(_mesh->regions()[r]);
		for (size_t i = 0; i < steps; i++) {
			step.step = i;
			step.currentTime = i;

			region->addSettings(i);
			files = store(_mesh->regions()[r]->name, step, region);
			region->clearData();
			_settings.push_back(std::make_pair(step, files));
		}
		delete region;
	}
}

void ResultStore::storeSubSolution(const Step &step, const std::vector<Solution*> &solution, const std::vector<std::pair<ElementType, Property> > &properties)
{
	if (_configuration.subsolution) {
		storeSolution(step, solution, properties);
	}
}

void ResultStore::storeSolution(const Step &step, const std::vector<Solution*> &solution, const std::vector<std::pair<ElementType, Property> > &properties)
{
	if (!_configuration.solution) {
		return;
	}

	prepare();

	_meshInfo->addSolution(solution);
	for (auto p = properties.begin(); p != properties.end(); ++p) {
		_meshInfo->addProperty(step, p->first, p->second);
	}
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
	if (!_configuration.FETI_data) {
		return;
	}
	if (_meshInfo == NULL) {
		if (_configuration.collected) {
			_meshInfo = new CollectedInfo(_mesh, _mode);
		} else {
			_meshInfo = new DistributedInfo(_mesh, _configuration.domain_shrink_ratio, _configuration.cluster_shrink_ratio, _mode);
		}
	}

	storeElementInfo(step);
	storeFixPoints(step);
	storeCorners(step);
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

	Region region(ElementType::NODES, fixPoints);

	MeshInfo *info = _meshInfo->deriveRegion(&region);
	std::vector<std::string> files = store("fix_points", step, info);
	delete info;
	_FETIdata.push_back(std::make_pair(step, files));
}

void ResultStore::storeCorners(const Step &step)
{
	std::vector<Element*> corners = _mesh->corners();

	Region region(ElementType::NODES, corners);

	MeshInfo *info = _meshInfo->deriveRegion(&region);
	std::vector<std::string> files = store("corners", step, info);
	delete info;
	_FETIdata.push_back(std::make_pair(step, files));
}

void ResultStore::storeDirichlet(const Step &step, const Instance &instance)
{
	if (!_meshInfo->distributed()) {
		// TODO: send dirichlet to process 0
		return;
	}

	for (size_t p = 0; p < instance.properties.size(); p++) {
		MeshInfo *info = _meshInfo->copyWithoutMesh();
		std::vector<double> *values = new std::vector<double>();
		Point point;

		for (size_t d = 0; d < instance.domains; d++) {
			for (size_t i = 0; i < instance.B1[d].I_row_indices.size() && instance.B1[d].I_row_indices[i] <= (eslocal)instance.block[Instance::CONSTRAINT::DIRICHLET]; i++) {
				const Element *e = _mesh->getDOFsElement(d, instance.B1[d].J_col_indices[i] - 1);
				if (e->DOFOffset(d, instance.B1[d].J_col_indices[i] - 1) != p) {
					continue;
				}
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
		info->_regions[0].data.elementDataDouble[ss.str()] = std::make_pair(1, values);
		std::vector<std::string> files = store("DIRICHLET_" + ss.str(), step, info);
		_FETIdata.push_back(std::make_pair(step, files));
		delete info;
	}

}

void ResultStore::storeLambdas(const Step &step, const Instance &instance)
{
	if (!_meshInfo->distributed()) {
		// it is pointless to store lambdas for collected result
		return;
	}

	std::vector<int> neighbours(environment->MPIsize);
	std::iota(neighbours.begin(), neighbours.end(), 0);

	for (size_t p = 0; p < instance.properties.size(); p++) {
		MeshInfo *info = _meshInfo->copyWithoutMesh();

		std::vector<double> duplicity, *values = new std::vector<double>();
		std::vector<std::pair<eslocal, eslocal> > lambdaMap;

		std::vector<std::vector<esglobal> > sLambdas(environment->MPIsize);
		std::vector<std::vector<Point> > sPoints(environment->MPIsize);
		std::vector<std::vector<esglobal> > rLambdas(environment->MPIsize);
		std::vector<std::vector<Point> > rPoints(environment->MPIsize);
		Point point;

		for (size_t d = 0; d < instance.domains; d++) {
			auto start = std::upper_bound(instance.B1[d].I_row_indices.begin(), instance.B1[d].I_row_indices.end(), instance.block[Instance::CONSTRAINT::DIRICHLET]);
			auto end = std::upper_bound(instance.B1[d].I_row_indices.begin(), instance.B1[d].I_row_indices.end(), instance.block[Instance::CONSTRAINT::DIRICHLET] + instance.block[Instance::CONSTRAINT::EQUALITY_CONSTRAINTS]);
			for (eslocal i = start - instance.B1[d].I_row_indices.begin(); i < end - instance.B1[d].I_row_indices.begin(); i++) {
				auto it = std::lower_bound(instance.B1clustersMap.begin(), instance.B1clustersMap.end(), instance.B1[d].I_row_indices[i] - 1, [&] (const std::vector<esglobal> &v, esglobal i) {
					return v[0] < i;
				});
				const Element *e = _mesh->getDOFsElement(d, instance.B1[d].J_col_indices[i] - 1);
				if (e->DOFOffset(d, instance.B1[d].J_col_indices[i] - 1) != p) {
					continue;
				}

				point = info->shrink(_mesh->coordinates()[e->node(0)], d);
				info->_regions[0].coordinates.insert(info->_regions[0].coordinates.end(), { point.x, point.y, point.z, 0, 0, 0 });
				lambdaMap.push_back(std::make_pair(instance.B1[d].I_row_indices[i], (eslocal)duplicity.size()));
				values->push_back(instance.B1[d].V_values[i]);
				duplicity.push_back(instance.B1duplicity[d][i]);

				if (it->size() == 3) {
					sLambdas[(*it)[2]].push_back(it->front());
					sPoints[(*it)[2]].push_back(point);
				}
			}
		}

		std::sort(lambdaMap.begin(), lambdaMap.end());


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

		for (size_t d = 0, offset = 0; d < instance.domains; d++) {
			auto start = std::upper_bound(instance.B1[d].I_row_indices.begin(), instance.B1[d].I_row_indices.end(), instance.block[Instance::CONSTRAINT::DIRICHLET]);
			auto end = std::upper_bound(instance.B1[d].I_row_indices.begin(), instance.B1[d].I_row_indices.end(), instance.block[Instance::CONSTRAINT::DIRICHLET] + instance.block[Instance::CONSTRAINT::EQUALITY_CONSTRAINTS]);
			for (eslocal i = start - instance.B1[d].I_row_indices.begin(); i < end - instance.B1[d].I_row_indices.begin(); i++) {
				auto it = std::lower_bound(instance.B1clustersMap.begin(), instance.B1clustersMap.end(), instance.B1[d].I_row_indices[i] - 1, [&] (const std::vector<esglobal> &v, esglobal i) {
					return v[0] < i;
				});
				const Element *e = _mesh->getDOFsElement(d, instance.B1[d].J_col_indices[i] - 1);
				if (e->DOFOffset(d, instance.B1[d].J_col_indices[i] - 1) != p) {
					continue;
				}

				Point a(info->_regions[0].coordinates[6 * offset + 0], info->_regions[0].coordinates[6 * offset + 1], info->_regions[0].coordinates[6 * offset + 2]);
				Point b;
				if (it->size() == 3) {
					size_t index = std::find(rLambdas[(*it)[2]].begin(), rLambdas[(*it)[2]].end(), it->front()) - rLambdas[(*it)[2]].begin();
					if (index == rLambdas[(*it)[2]].size() || rLambdas[(*it)[2]][index] != it->front()) {
						ESINFO(ERROR) << "Different Lambdas on neighbour clusters.";
					}
					b = rPoints[(*it)[2]][index];
				} else {
					auto lit = std::lower_bound(lambdaMap.begin(), lambdaMap.end(), std::make_pair(instance.B1[d].I_row_indices[i], (eslocal)-1));
					if (lit == lambdaMap.end()) {
						ESINFO(ERROR) << "Incorrect B1 lambdas.";
					}
					if (lit->second == (eslocal)offset) {
						lit++;
					}
					b.x = info->_regions[0].coordinates[6 * lit->second + 0];
					b.y = info->_regions[0].coordinates[6 * lit->second + 1];
					b.z = info->_regions[0].coordinates[6 * lit->second + 2];
				}
				Point length = b - a;
				Point ax = a + length * duplicity[offset];
				info->_regions[0].coordinates[6 * offset + 3] = ax.x;
				info->_regions[0].coordinates[6 * offset + 4] = ax.y;
				info->_regions[0].coordinates[6 * offset + 5] = ax.z;
				offset++;
			}
		}


		size_t size = info->_regions[0].coordinates.size() / 6;
		info->_regions[0].elementsTypes.resize(size, Line2VTKCode);
		info->_regions[0].elementsNodes.reserve(size);
		info->_regions[0].elements.reserve(2 * size);

		for (size_t i = 0; i < size; i++) {
			info->_regions[0].elementsNodes.push_back(2 * i + 2);
			info->_regions[0].elements.push_back(2 * i);
			info->_regions[0].elements.push_back(2 * i + 1);
		}

		info->_regions[0].data.elementDataDouble["B1values"] = std::make_pair(1, values);
		std::vector<std::string> files = store("EQUALITY_CONSTRAINTS", step, info);
		_FETIdata.push_back(std::make_pair(step, files));
		delete info;
	}
}






