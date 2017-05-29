
#include "gridtower.h"

#include "../../../configuration/environment.h"
#include "../../../configuration/input/inputgeneratorgridtower.h"

#include "../../../mesh/structures/mesh.h"
#include "../../../mesh/structures/elementtypes.h"
#include "../primitives/block.h"

#include "mpi.h"

using namespace espreso::input;

void GridTower::load(const GridTowerConfiguration &configuration, Mesh &mesh, size_t index, size_t size)
{
	GridTower gridtower(configuration, mesh, index, size);
	gridtower.fill();
}

GridTower::GridTower(const GridTowerConfiguration &configuration, Mesh &mesh, size_t index, size_t size)
: Loader(mesh), _gridTower(configuration), _grid(NULL), _gridPointsIDOffset(0), _index(index), _size(size)
{
	size_t gridIndex = 0, firstCluster = 0, lastCluster = 0;
	Triple<double> gridPointsOffset(0, 0, 0);
	for (auto it = _gridTower.grids.begin(); it != _gridTower.grids.end(); ++it, gridIndex++) {
		lastCluster += it->second->clusters_x * it->second->clusters_y * it->second->clusters_z;
		if (firstCluster <= _index && _index < lastCluster) {
			_grid = new Grid(*it->second, mesh, _index - firstCluster, lastCluster - firstCluster);
			_clusterIndexBegin = firstCluster;
			_clusterIndexEnd = lastCluster;
			_gridIndex = gridIndex;
			_gridPointsOffset = gridPointsOffset;
		}
		firstCluster = lastCluster;
		gridPointsOffset += Triple<double>(it->second->length_x, it->second->length_y, it->second->length_z);
	}

	if (lastCluster != size) {
		ESINFO(GLOBAL_ERROR) << "Incorrect number of MPI processes (" << size << "). Should be " << lastCluster;
	}

	size_t gridPointsCounter = _clusterIndexBegin == index ? _grid->pointCount() : 0;

	std::vector<size_t> counters(size);
	MPI_Allgather(&gridPointsCounter, sizeof(size_t), MPI_BYTE, counters.data(), sizeof(size_t), MPI_BYTE, environment->MPICommunicator);
	for (size_t i = 0; i < _clusterIndexBegin; i++) {
		_gridPointsIDOffset += counters[i];
	}
	_grid->bodyIndex(_gridIndex);
}

GridTower::~GridTower()
{
	delete _grid;
}

void GridTower::points(Coordinates &coordinates)
{
	Triple<double> end = _grid->_block->block.end, start = _grid->_block->block.start;
	switch (_gridTower.direction) {
	case GridTowerConfiguration::DIRECTION::X:
		_grid->_block->block.start.x += _gridPointsOffset.x;
		_grid->_block->block.end.x += _gridPointsOffset.x;
		break;
	case GridTowerConfiguration::DIRECTION::Y:
		_grid->_block->block.start.y += _gridPointsOffset.y;
		_grid->_block->block.end.y += _gridPointsOffset.y;
		break;
	case GridTowerConfiguration::DIRECTION::Z:
		_grid->_block->block.start.z += _gridPointsOffset.z;
		_grid->_block->block.end.z += _gridPointsOffset.z;
		break;
	}

	_grid->points(coordinates, _gridPointsIDOffset);

	_grid->_block->block.end = end;
	_grid->_block->block.start = start;
}

void GridTower::elements(std::vector<size_t> &bodies, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges)
{
	_grid->elements(bodies, elements, faces, edges);
	bodies.resize(_gridTower.grids.size() + 1);
	for (size_t i = 0; i < bodies.size(); i++) {
		if (i <= _gridIndex) {
			bodies[i] = 0;
		} else {
			bodies[i] = elements.size();
		}
	}
}

bool GridTower::partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners)
{
	return _grid->partitiate(nodes, partsPtrs, fixPoints, corners);
}

void GridTower::neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours, const std::vector<Element*> &faces, const std::vector<Element*> &edges)
{
	auto middle = _gridTower.grids.find(_gridIndex);
	auto lower  = _gridTower.grids.end();
	auto upper  = _gridTower.grids.end();

	switch (_gridTower.direction) {
	case GridTowerConfiguration::DIRECTION::X:
		if (_grid->_clusterOffset.x == 0) {
			lower  = _gridTower.grids.find(_gridIndex - 1);
		}
		if (_grid->_clusterOffset.x == _grid->_settings.clusters.x - 1) {
			upper  = _gridTower.grids.find(_gridIndex + 1);
		}
		break;
	case GridTowerConfiguration::DIRECTION::Y:
		if (_grid->_clusterOffset.y == 0) {
			lower  = _gridTower.grids.find(_gridIndex - 1);
		}
		if (_grid->_clusterOffset.y == _grid->_settings.clusters.y - 1) {
			upper  = _gridTower.grids.find(_gridIndex + 1);
		}
		break;
	case GridTowerConfiguration::DIRECTION::Z:
		if (_grid->_clusterOffset.z == 0) {
			lower  = _gridTower.grids.find(_gridIndex - 1);
		}
		if (_grid->_clusterOffset.z == _grid->_settings.clusters.z - 1) {
			upper  = _gridTower.grids.find(_gridIndex + 1);
		}
		break;
	}

	if (lower != _gridTower.grids.end()) {
		if (lower->second->blocks_z != middle->second->blocks_z || lower->second->blocks_y != middle->second->blocks_y || lower->second->blocks_x != middle->second->blocks_x) {
			ESINFO(GLOBAL_ERROR) << "Implement neighbors for grid tower with different number of blocks.";
		}
		if (lower->second->clusters_z != middle->second->clusters_z || lower->second->clusters_y != middle->second->clusters_y || lower->second->clusters_x != middle->second->clusters_x) {
			ESINFO(GLOBAL_ERROR) << "Implement neighbors for grid tower with different number of clusters.";
		}
		if (lower->second->domains_z != middle->second->domains_z || lower->second->domains_y != middle->second->domains_y || lower->second->domains_x != middle->second->domains_x) {
			ESINFO(GLOBAL_ERROR) << "Implement neighbors for grid tower with different number of domains.";
		}
		if (lower->second->elements_z != middle->second->elements_z || lower->second->elements_y != middle->second->elements_y || lower->second->elements_x != middle->second->elements_x) {
			ESINFO(GLOBAL_ERROR) << "Implement neighbors for grid tower with different number of elements.";
		}
	}
	if (upper != _gridTower.grids.end()) {
		if (upper->second->blocks_z != middle->second->blocks_z || upper->second->blocks_y != middle->second->blocks_y || upper->second->blocks_x != middle->second->blocks_x) {
			ESINFO(GLOBAL_ERROR) << "Implement neighbors for grid tower with different number of blocks.";
		}
		if (upper->second->clusters_z != middle->second->clusters_z || upper->second->clusters_y != middle->second->clusters_y || upper->second->clusters_x != middle->second->clusters_x) {
			ESINFO(GLOBAL_ERROR) << "Implement neighbors for grid tower with different number of clusters.";
		}
		if (upper->second->domains_z != middle->second->domains_z || upper->second->domains_y != middle->second->domains_y || upper->second->domains_x != middle->second->domains_x) {
			ESINFO(GLOBAL_ERROR) << "Implement neighbors for grid tower with different number of domains.";
		}
		if (upper->second->elements_z != middle->second->elements_z || upper->second->elements_y != middle->second->elements_y || upper->second->elements_x != middle->second->elements_x) {
			ESINFO(GLOBAL_ERROR) << "Implement neighbors for grid tower with different number of elements.";
		}
	}

	std::vector<int> map(27);

	Triple<int> offset, coffset(_grid->_clusterOffset), count = _grid->_settings.blocks * _grid->_settings.clusters, size = count.toSize();
	size_t index = 0;
	for (offset.z = -1; offset.z <= 1; offset.z++) {
		for (offset.y = -1; offset.y <= 1; offset.y++) {
			for (offset.x = -1; offset.x <= 1; offset.x++, index++) {
				if ((coffset + offset) < 0 || (count - (coffset + offset)) < 1) {
					map[index] = -1;
				} else {
					map[index] = _grid->_cMap[((coffset + offset) * size).sum()] + _clusterIndexBegin;
					if (map[index] != environment->MPIrank) {
						neighbours.push_back(map[index]);
					}
				}
			}
		}
	}


	size_t cx = middle->second->blocks_x * middle->second->clusters_x;
	size_t cy = middle->second->blocks_y * middle->second->clusters_y;
	size_t cz = middle->second->blocks_z * middle->second->clusters_z;

	for (offset.z = -1, index = 0; offset.z <= 1; offset.z++) {
		for (offset.y = -1; offset.y <= 1; offset.y++) {
			for (offset.x = -1; offset.x <= 1; offset.x++, index++) {
				switch (_gridTower.direction) {
				case GridTowerConfiguration::DIRECTION::X:
					if (offset.x == 1 && map[index - 1] != -1 && upper != _gridTower.grids.end()) {
						map[index] = map[index - 1] + cx * cy * (cz - 1) + cx * (cy - 1) + 1;
						neighbours.push_back(map[index]);
					}
					if (offset.x == -1 && map[index + 1] != -1 && lower != _gridTower.grids.end()) {
						map[index] = map[index + 1] - cx * cy * (cz - 1) - cx * (cy - 1) - 1;
						neighbours.push_back(map[index]);
					}
					break;
				case GridTowerConfiguration::DIRECTION::Y:
					if (offset.y == 1 && map[index - 3] != -1 && upper != _gridTower.grids.end()) {
						map[index] = map[index - 3] + cx * cy * (cz - 1) + cx;
						neighbours.push_back(map[index]);
					}
					if (offset.y == -1 && map[index + 3] != -1 && lower != _gridTower.grids.end()) {
						map[index] = map[index + 3] - cx * cy * (cz - 1) - cx;
						neighbours.push_back(map[index]);
					}
					break;
				case GridTowerConfiguration::DIRECTION::Z:
					if (offset.z == 1 && map[index - 9] != -1 && upper != _gridTower.grids.end()) {
						map[index] = map[index - 9] + cx * cy;
						neighbours.push_back(map[index]);
					}
					if (offset.z == -1 && map[index + 9] != -1 && lower != _gridTower.grids.end()) {
						map[index] = map[index + 9] - cx * cy;
						neighbours.push_back(map[index]);
					}
					break;
				}
			}
		}
	}

	_grid->_block->boundaries(nodes, map);
	std::sort(neighbours.begin(), neighbours.end());
	mesh.synchronizeGlobalIndices();
}

void GridTower::regions(
		std::vector<Evaluator*> &evaluators,
		std::vector<Region*> &regions,
		std::vector<Element*> &elements,
		std::vector<Element*> &faces,
		std::vector<Element*> &edges,
		std::vector<Element*> &nodes)
{
	_grid->regions(evaluators, regions, elements, faces, edges, nodes);

	for (auto it = _gridTower.grids.begin(); it != _gridTower.grids.end(); ++it) {
		for (auto r = it->second->nodes.begin(); r != it->second->nodes.end(); ++r) {
			auto it = std::find_if(regions.begin(), regions.end(), [&] (const Region *region) { return region->name.compare(r->first) == 0; });
			if (it == mesh.regions().end()) {
				regions.push_back(new Region(ElementType::NODES));
				regions.back()->name = r->first;
			}
		}
	}
	mesh.synchronizeRegionOrder();
}



