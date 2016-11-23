
#include "grid.h"

#include "../../../config/description.h"

using namespace espreso::input;

Grid::Grid(Mesh &mesh, size_t index, size_t size)
: Loader(mesh), _index(index), _size(size)
{
	const GridGenerator &grid = configuration.generator.grid;

	_generator.blocks     = Triple<size_t>(grid.blocks_x, grid.blocks_y, grid.blocks_z);
	_generator.clusters = Triple<size_t>(grid.clusters_x, grid.clusters_y, grid.clusters_z);
	_generator.nonempty.resize((_generator.blocks * _generator.clusters).mul(), true);

	_generator.domains  = Triple<size_t>(grid.domains_x, grid.domains_y, grid.domains_z);
	_generator.elements = Triple<size_t>(grid.elements_x, grid.elements_y, grid.elements_z);

	_generator.start = Triple<double>(grid.start_x, grid.start_y, grid.start_z);
	_generator.end   = Triple<double>(grid.start_x + grid.length_x, grid.start_y + grid.length_y, grid.start_z + grid.length_z);

	for (auto it = grid.blocks.values.begin(); it != grid.blocks.values.end(); ++it) {
		if (it->first >= _generator.nonempty.size()) {
			ESINFO(GLOBAL_ERROR) << "Block index is out of range.";
		}
		_generator.nonempty[it->first] = it->second;
	}

	Triple<size_t> clusters = _generator.blocks * _generator.clusters;

	size_t cluster = 0;
	Triple<size_t> offset;
	int clusterIndex = 0;
	for (offset.z = 0; offset.z < clusters.z; offset.z++) {
		for (offset.y = 0; offset.y < clusters.y; offset.y++) {
			for (offset.x = 0; offset.x < clusters.x; offset.x++, clusterIndex++) {

				if (_generator.nonempty[cluster]) {
					_cMap.push_back(cluster);
				} else {
					_cMap.push_back(-1);
					continue;
				}

				if (cluster++ == index) {
					_clusterOffset = offset;
					_clusterIndex = clusterIndex;
					BlockSetting block;
					block.domains  = _generator.domains;
					block.elements = _generator.elements;
					block.start = _generator.start + (_generator.end - _generator.start) / _generator.clusters * offset;
					block.end   = _generator.start + (_generator.end - _generator.start) / _generator.clusters * (offset + 1);
					_block = new Block<Hexahedron8>(mesh, block);
					_subnodes = Hexahedron8::subnodes;
				}
			}
		}
	}
	if (cluster != size) {
		ESINFO(GLOBAL_ERROR) << "Incorrect number of MPI processes (" << size << "). Should be " << cluster;
	}
}

void Grid::points(Coordinates &coordinates)
{
	_block->points(coordinates._points);

	Triple<size_t> cnodes = _generator.domains * _generator.elements * (Triple<size_t>(_subnodes) + 1);
	Triple<size_t> gnodes = _generator.blocks * _generator.clusters * cnodes;
	Triple<size_t> coffset = _clusterOffset * cnodes;
	Triple<size_t> size = (gnodes + 1).toSize();


	Triple<size_t> offset;
	for (offset.z = 0; offset.z <= cnodes.z; offset.z++) {
		for (offset.y = 0; offset.y <= cnodes.y; offset.y++) {
			for (offset.x = 0; offset.x <= cnodes.x; offset.x++) {
				coordinates._globalIndex.push_back(((coffset + offset) * size).sum());
				coordinates._globalMapping.push_back(std::make_pair(coordinates._globalIndex.back(), (esglobal)coordinates._globalMapping.size()));
			}
		}
	}
}

void Grid::elements(std::vector<Element*> &elements)
{
	_block->elements(elements);
}

void Grid::materials(std::vector<Material> &materials)
{
	materials.push_back(Material(mesh.coordinates()));
}

void Grid::clusterBoundaries(std::vector<Element*> &nodes, std::vector<int> &neighbours)
{
	std::vector<int> map(27);

	Triple<int> offset, coffset(_clusterOffset), count = _generator.blocks * _generator.clusters, size = count.toSize();
	size_t index = 0;
	for (offset.z = -1; offset.z <= 1; offset.z++) {
		for (offset.y = -1; offset.y <= 1; offset.y++) {
			for (offset.x = -1; offset.x <= 1; offset.x++, index++) {
				if ((coffset + offset) < 0 || (count - (coffset + offset)) < 1) {
					map[index] = -1;
				} else {
					map[index] = _cMap[((coffset + offset) * size).sum()];
					if (map[index] != config::env::MPIrank) {
						neighbours.push_back(map[index]);
					}
				}
			}
		}
	}

	_block->boundaries(nodes, map);
	std::sort(neighbours.begin(), neighbours.end());
}

void Grid::settings(
		std::vector<Evaluator*> &evaluators,
		std::vector<Region> &regions,
		std::vector<Element*> &elements,
		std::vector<Element*> &faces,
		std::vector<Element*> &edges,
		std::vector<Element*> &nodes)
{
	for (auto it = configuration.generator.grid.nodes.values.begin(); it != configuration.generator.grid.nodes.values.end(); ++it) {
		regions.push_back(Region());
		regions.back().name = it->first;
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			regions.back().elements = nodes;
		} else {
			BlockBorder border(it->second);
			_block->region(nodes, regions.back(), border, 0);
		}
	}
	for (auto it = configuration.generator.grid.edges.values.begin(); it != configuration.generator.grid.edges.values.end(); ++it) {
		regions.push_back(Region());
		regions.back().name = it->first;
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			ESINFO(GLOBAL_ERROR) << "Implement region of all edges.";
		} else {
			BlockBorder border(it->second);
			_block->region(nodes, regions.back(), border, 1);
		}
	}
	for (auto it = configuration.generator.grid.faces.values.begin(); it != configuration.generator.grid.faces.values.end(); ++it) {
		regions.push_back(Region());
		regions.back().name = it->first;
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			ESINFO(GLOBAL_ERROR) << "Implement region of all faces.";
		} else {
			BlockBorder border(it->second);
			_block->region(nodes, regions.back(), border, 2);
		}
	}

	for (auto it = configuration.generator.grid.elements.values.begin(); it != configuration.generator.grid.elements.values.end(); ++it) {
		regions.push_back(Region());
		regions.back().name = it->first;
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			regions.back().elements = elements;
		} else {
			BlockBorder border(it->second);
			_block->region(nodes, regions.back(), border, 3);
		}
	}
}
