
#include "grid.h"
#include "../primitives/block.h"
#include "../elements/elements.h"

#include "../../../config/inputgeneratorgrid.h"

using namespace espreso::input;

GridSettings::GridSettings()
: etype(ELEMENT_TYPE::HEXA8),
  start(0, 0, 0), end(1, 1, 1),
  blocks(1, 1, 1), clusters(1, 1, 1), domains(2, 2, 2), elements(5, 5, 5),
  projection(Expression("x", { "x", "y", "z" }), Expression("y", { "x", "y", "z" }), Expression("z", { "x", "y", "z" })),
  rotation(Expression("0", { "x", "y", "z" }), Expression("0", { "x", "y", "z" }), Expression("0", { "x", "y", "z" })),
  nonempty(1, true), uniformDecomposition(true) {}

GridSettings::GridSettings(const GridConfiguration &configuration)
: projection(Expression(configuration.projection_x, { "x", "y", "z" }), Expression(configuration.projection_y, { "x", "y", "z" }), Expression(configuration.projection_z, { "x", "y", "z" })),
  rotation(Expression(configuration.rotation_x, { "x", "y", "z" }), Expression(configuration.rotation_y, { "x", "y", "z" }), Expression(configuration.rotation_z, { "x", "y", "z" }))
{
	etype = configuration.element_type;

	blocks     = Triple<size_t>(configuration.blocks_x, configuration.blocks_y, configuration.blocks_z);
	clusters = Triple<size_t>(configuration.clusters_x, configuration.clusters_y, configuration.clusters_z);
	nonempty.resize((blocks * clusters).mul(), true);

	domains  = Triple<size_t>(configuration.domains_x, configuration.domains_y, configuration.domains_z);
	elements = Triple<size_t>(configuration.elements_x, configuration.elements_y, configuration.elements_z);

	start = Triple<double>(configuration.start_x, configuration.start_y, configuration.start_z);
	end   = Triple<double>(configuration.start_x + configuration.length_x, configuration.start_y + configuration.length_y, configuration.start_z + configuration.length_z);

	for (auto it = configuration.blocks.values.begin(); it != configuration.blocks.values.end(); ++it) {
		if (it->first >= nonempty.size()) {
			ESINFO(GLOBAL_ERROR) << "Block index is out of range.";
		}
		nonempty[it->first] = it->second;
	}

	uniformDecomposition = configuration.uniform_decomposition;
}

void Grid::load(const GridConfiguration &configuration, Mesh &mesh, size_t index, size_t size)
{
	ESINFO(OVERVIEW) << "Generate grid";
	Grid grid(configuration, mesh, index, size);
	grid.fill();
}

Grid::Grid(const GridConfiguration &configuration, Mesh &mesh, size_t index, size_t size)
: Loader(mesh), _grid(configuration), _settings(configuration), _index(index), _size(size)
{
	Triple<size_t> clusters = _settings.blocks * _settings.clusters;

	size_t cluster = 0;
	Triple<size_t> offset;
	int clusterIndex = 0;
	for (offset.z = 0; offset.z < clusters.z; offset.z++) {
		for (offset.y = 0; offset.y < clusters.y; offset.y++) {
			for (offset.x = 0; offset.x < clusters.x; offset.x++, clusterIndex++) {

				if (_settings.nonempty[cluster]) {
					_cMap.push_back(cluster);
				} else {
					_cMap.push_back(-1);
					continue;
				}

				if (cluster++ == index) {
					_clusterOffset = offset;
					_clusterIndex = clusterIndex;
					BlockSetting block;
					block.domains  = _settings.domains;
					block.elements = _settings.elements;
					block.start = _settings.start + (_settings.end - _settings.start) / _settings.clusters * offset;
					block.end   = _settings.start + (_settings.end - _settings.start) / _settings.clusters * (offset + 1);
					block.projection = _settings.projection;
					block.rotation = _settings.rotation;
					switch (_settings.etype) {
					case ELEMENT_TYPE::HEXA8:
						_block = new Block<Hexahedron8>(mesh, block);
						_subnodes = Hexahedron8::subnodes;
						break;
					case ELEMENT_TYPE::HEXA20:
						_block = new Block<Hexahedron20>(mesh, block);
						_subnodes = Hexahedron20::subnodes;
						break;
					case ELEMENT_TYPE::TETRA4:
						_block = new Block<Tetrahedron4>(mesh, block);
						_subnodes = Tetrahedron4::subnodes;
						break;
					case ELEMENT_TYPE::TETRA10:
						_block = new Block<Tetrahedron10>(mesh, block);
						_subnodes = Tetrahedron10::subnodes;
						break;
					case ELEMENT_TYPE::PRISMA6:
						_block = new Block<Prisma6>(mesh, block);
						_subnodes = Prisma6::subnodes;
						break;
					case ELEMENT_TYPE::PRISMA15:
						_block = new Block<Prisma15>(mesh, block);
						_subnodes = Prisma15::subnodes;
						break;
					case ELEMENT_TYPE::PYRAMID5:
						_block = new Block<Pyramid5>(mesh, block);
						_subnodes = Pyramid5::subnodes;
						break;
					case ELEMENT_TYPE::PYRAMID13:
						_block = new Block<Pyramid13>(mesh, block);
						_subnodes = Pyramid13::subnodes;
						break;

					case ELEMENT_TYPE::SQUARE4:
						_block = new Block<Square4>(mesh, block);
						_subnodes = Square4::subnodes;
						break;
					case ELEMENT_TYPE::SQUARE8:
						_block = new Block<Square8>(mesh, block);
						_subnodes = Square8::subnodes;
						break;
					case ELEMENT_TYPE::TRIANGLE3:
						_block = new Block<Triangle3>(mesh, block);
						_subnodes = Triangle3::subnodes;
						break;
					case ELEMENT_TYPE::TRIANGLE6:
						_block = new Block<Triangle6>(mesh, block);
						_subnodes = Triangle6::subnodes;
						break;
					}
				}
			}
		}
	}
	if (cluster != size) {
		ESINFO(GLOBAL_ERROR) << "Incorrect number of MPI processes (" << size << "). Should be " << cluster;
	}
}

Grid::~Grid()
{
	delete _block;
}

void Grid::points(Coordinates &coordinates)
{
	_block->points(coordinates._points);

	Triple<size_t> cnodes = _settings.domains * _settings.elements * (Triple<size_t>(_subnodes) - 1);
	Triple<size_t> gnodes = _settings.blocks * _settings.clusters * cnodes;
	Triple<size_t> coffset = _clusterOffset * cnodes;
	Triple<size_t> size = (gnodes + 1).toSize();


	Triple<size_t> offset;
	for (offset.z = 0; offset.z <= cnodes.z; offset.z++) {
		for (offset.y = 0; offset.y <= cnodes.y; offset.y++) {
			for (offset.x = 0; offset.x <= cnodes.x; offset.x++) {
				coordinates._globalIndex.push_back(((coffset + offset) * size).sum());
				coordinates._globalMapping.push_back(G2L(coordinates._globalIndex.back(), (esglobal)coordinates._globalMapping.size()));
			}
		}
	}
}

void Grid::elements(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges)
{
	_block->elements(elements);
}

bool Grid::partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners)
{
	if (_settings.uniformDecomposition) {
		_block->uniformPartition(partsPtrs, _settings.domains.mul());
		_block->uniformFixPoints(nodes, fixPoints);
		_block->uniformCorners(nodes, corners, 1, true, true, true);
		return true;
	} else {
		mesh.partitiate(_settings.domains.mul());
		return false;
	}
}

void Grid::neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours, const std::vector<Element*> &faces, const std::vector<Element*> &edges)
{
	std::vector<int> map(27);

	Triple<int> offset, coffset(_clusterOffset), count = _settings.blocks * _settings.clusters, size = count.toSize();
	size_t index = 0;
	for (offset.z = -1; offset.z <= 1; offset.z++) {
		for (offset.y = -1; offset.y <= 1; offset.y++) {
			for (offset.x = -1; offset.x <= 1; offset.x++, index++) {
				if ((coffset + offset) < 0 || (count - (coffset + offset)) < 1) {
					map[index] = -1;
				} else {
					map[index] = _cMap[((coffset + offset) * size).sum()];
					if (map[index] != environment->MPIrank) {
						neighbours.push_back(map[index]);
					}
				}
			}
		}
	}

	_block->boundaries(nodes, map);
	std::sort(neighbours.begin(), neighbours.end());
}

void Grid::regions(
		std::vector<Evaluator*> &evaluators,
		std::vector<Region*> &regions,
		std::vector<Element*> &elements,
		std::vector<Element*> &faces,
		std::vector<Element*> &edges,
		std::vector<Element*> &nodes)
{
	for (auto it = _grid.nodes.values.begin(); it != _grid.nodes.values.end(); ++it) {
		regions.push_back(new Region());
		regions.back()->name = it->first;
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			regions.back()->elements = nodes;
		} else {
			BlockBorder border(it->second);
			_block->region(nodes, regions.back(), border, 0);
		}
	}
	for (auto it = _grid.edges.values.begin(); it != _grid.edges.values.end(); ++it) {
		regions.push_back(new Region());
		regions.back()->name = it->first;
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			ESINFO(GLOBAL_ERROR) << "Implement region of all edges.";
		} else {
			BlockBorder border(it->second);
			_block->region(nodes, regions.back(), border, 1);
			edges.insert(edges.end(), regions.back()->elements.begin(), regions.back()->elements.end());
		}
	}
	for (auto it = _grid.faces.values.begin(); it != _grid.faces.values.end(); ++it) {
		regions.push_back(new Region());
		regions.back()->name = it->first;
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			ESINFO(GLOBAL_ERROR) << "Implement region of all faces.";
		} else {
			BlockBorder border(it->second);
			_block->region(nodes, regions.back(), border, 2);
			faces.insert(faces.end(), regions.back()->elements.begin(), regions.back()->elements.end());
		}
	}

	for (auto it = _grid.elements.values.begin(); it != _grid.elements.values.end(); ++it) {
		regions.push_back(new Region());
		regions.back()->name = it->first;
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			regions.back()->elements = elements;
		} else {
			BlockBorder border(it->second);
			_block->region(elements, regions.back(), border, 3);
		}
	}
}
