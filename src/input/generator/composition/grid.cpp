
#include "grid.h"

using namespace espreso::input;

Grid::Grid(Mesh &mesh, size_t index, size_t size)
: Loader(mesh), _index(index), _size(size)
{
	_generator.grid     = Triple<size_t>(1, 1, 1);
	_generator.clusters = Triple<size_t>(1, 1, 1);
	_generator.nonempty.resize((_generator.grid * _generator.clusters).mul(), true);

	_generator.domains  = Triple<size_t>(2, 2, 2);
	_generator.elements = Triple<size_t>(5, 5, 5);

	_generator.start = Triple<double>(0, 0, 0);
	_generator.end   = Triple<double>(30, 30, 30);

	Triple<size_t> clusters = _generator.grid * _generator.clusters;

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
	Triple<size_t> gnodes = _generator.grid * _generator.clusters * cnodes;
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

	Triple<int> offset, coffset(_clusterOffset), count = _generator.grid * _generator.clusters, size = count.toSize();
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
	BlockBorder bottom; bottom.end.x = 30; bottom.end.y = 30;
	BlockBorder left; left.end.x = 30; left.end.z = 30;
	BlockBorder back; back.end.y = 30; back.end.z = 30;

	BlockBorder right; right.end.x = 30; right.end.z = 30; right.end.z = right.start.z = 30;
	BlockBorder top; top.end.x = 30; top.end.z = 30; top.end.z = top.start.z = 30;
	BlockBorder front; front.end.x = 30; front.end.z = 30; front.end.z = front.start.z = 30;

	Region region1;
	_block->region(nodes, region1, bottom, 0);
	evaluators.push_back(new ConstEvaluator(0));
	for (size_t i = 0; i < region1.elements.size(); i++) {
		region1.elements[i]->addSettings(Property::DISPLACEMENT_Z, evaluators.back());
	}

	Region region2;
	_block->region(nodes, region2, left, 0);
	for (size_t i = 0; i < region2.elements.size(); i++) {
		region2.elements[i]->addSettings(Property::DISPLACEMENT_Y, evaluators.back());
	}

	Region region3;
	_block->region(nodes, region3, back, 0);
	for (size_t i = 0; i < region3.elements.size(); i++) {
		region3.elements[i]->addSettings(Property::DISPLACEMENT_X, evaluators.back());
	}
}
