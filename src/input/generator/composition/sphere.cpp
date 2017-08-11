

#include "sphere.h"

#include "../../../configuration/input/inputgeneratorsphere.h"
#include "../primitives/block.h"
#include "../generator.h"
#include "../../../mesh/structures/region.h"
#include "../../../mesh/structures/mesh.h"
#include "../../../mesh/structures/elementtypes.h"
#include "../../../mesh/structures/coordinates.h"

#include "../elements/3D/hexahedron20.h"
#include "../elements/3D/hexahedron8.h"
#include "../elements/3D/prisma15.h"
#include "../elements/3D/prisma6.h"
#include "../elements/3D/pyramid13.h"
#include "../elements/3D/pyramid5.h"
#include "../elements/3D/tetrahedron10.h"
#include "../elements/3D/tetrahedron4.h"

using namespace espreso::input;

SphereSettings::SphereSettings()
: etype(ELEMENT_TYPE::HEXA8),
  innerRadius(5), outerRadius(10),
  clusters(1), layers(1), domains(2, 2, 2), elements(5, 5, 5),
  uniformDecomposition(true) {}

SphereSettings::SphereSettings(const SphereConfiguration &configuration)
{
	etype = configuration.element_type;

	clusters = configuration.clusters;
	layers = configuration.layers;

	domains  = Triple<size_t>(configuration.domains_x, configuration.domains_y, configuration.domains_z);
	elements = Triple<size_t>(configuration.elements_x, configuration.elements_y, configuration.elements_z);

	innerRadius = configuration.inner_radius;
	outerRadius = configuration.outer_radius;

	uniformDecomposition = configuration.uniform_decomposition;
}

void Sphere::load(const SphereConfiguration &configuration, Mesh &mesh, size_t index, size_t size)
{
	Sphere sphere(configuration, mesh, index, size);
	sphere.fill();
}

Sphere::Sphere(const SphereConfiguration &configuration, Mesh &mesh, size_t index, size_t size)
: Loader(mesh), _sphere(configuration), _settings(configuration), _index(index), _size(size)
{
	if (size % 6 != 0) {
		ESINFO(GLOBAL_ERROR) << "Number of MPI process should be 6 x clusters for sphere generator.";
	}

	BlockSetting block;
	block.domains  = _settings.domains;
	block.elements = _settings.elements;

	_row = ((_index % (_settings.clusters * _settings.clusters)) % _settings.clusters);
	_col = ((_index % (_settings.clusters * _settings.clusters)) / _settings.clusters);
	_layer = _index / (6 * _settings.clusters * _settings.clusters);

	block.start = Triple<esglobal>( _row      / _settings.clusters / Generator::precision,  _col      / _settings.clusters / Generator::precision,  _layer      / _settings.layers / Generator::precision);
	block.end   = Triple<esglobal>((_row + 1) / _settings.clusters / Generator::precision, (_col + 1) / _settings.clusters / Generator::precision, (_layer + 1) / _settings.layers / Generator::precision);

	switch ((index / (_settings.clusters * _settings.clusters)) % 6) {
	case 0:
		_side = SIDE::UP;
		block.projection = Triple<Expression>(
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(x * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(y * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4))"                           , { "x", "y", "z" }));
		break;
	case 1:
		_side = SIDE::FRONT;
		block.projection = Triple<Expression>(
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4))"                           , { "x", "y", "z" }),
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(y * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(x * pi / 2 - pi / 4)", { "x", "y", "z" }));;
		break;
	case 2:
		_side = SIDE::DOWN;
		block.projection = Triple<Expression>(
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(x * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(y * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4))"                           , { "x", "y", "z" }));
		break;
	case 3:
		_side = SIDE::BACK;
		block.projection = Triple<Expression>(
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4))"                           , { "x", "y", "z" }),
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(y * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(x * pi / 2 - pi / 4)", { "x", "y", "z" }));;
		break;
	case 4:
		_side = SIDE::LEFT;
		block.projection = Triple<Expression>(
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(y * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4))"                           , { "x", "y", "z" }),
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(x * pi / 2 - pi / 4)", { "x", "y", "z" }));
		break;
	case 5:
		_side = SIDE::RIGHT;
		block.projection = Triple<Expression>(
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(y * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4))"                           , { "x", "y", "z" }),
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(x * pi / 2 - pi / 4)", { "x", "y", "z" }));
		break;
	}

	std::string element;

	switch (_settings.etype) {
	case ELEMENT_TYPE::HEXA8:
		_block = new Block<Hexahedron8>(block);
		_subnodes = Hexahedron8::subnodes;
		element = "hexahedrons";
		break;
	case ELEMENT_TYPE::HEXA20:
		_block = new Block<Hexahedron20>(block);
		_subnodes = Hexahedron20::subnodes;
		element = "hexahedrons with midnodes";
		break;
	case ELEMENT_TYPE::TETRA4:
		_block = new Block<Tetrahedron4>(block);
		_subnodes = Tetrahedron4::subnodes;
		element = "tetrahedrons";
		break;
	case ELEMENT_TYPE::TETRA10:
		_block = new Block<Tetrahedron10>(block);
		_subnodes = Tetrahedron10::subnodes;
		element = "tetrahedrons with midnodes";
		break;
	case ELEMENT_TYPE::PRISMA6:
		_block = new Block<Prisma6>(block);
		_subnodes = Prisma6::subnodes;
		element = "prismas";
		break;
	case ELEMENT_TYPE::PRISMA15:
		_block = new Block<Prisma15>(block);
		_subnodes = Prisma15::subnodes;
		element = "prismas with midnodes";
		break;
	case ELEMENT_TYPE::PYRAMID5:
		_block = new Block<Pyramid5>(block);
		_subnodes = Pyramid5::subnodes;
		element = "pyramids";
		break;
	case ELEMENT_TYPE::PYRAMID13:
		_block = new Block<Pyramid13>(block);
		_subnodes = Pyramid13::subnodes;
		element = "pyramids with midnodes";
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unsupported element type";
	}

	ESINFO(OVERVIEW) << "Generate sphere " << element;
}

Sphere::~Sphere()
{
	delete _block;
}

void Sphere::points(Coordinates &coordinates)
{
	_block->points(coordinates._points);

	for (size_t i = 0; i < coordinates._points.size(); i++) {
		coordinates._globalIndex.push_back(i + _index * coordinates._points.size());
		coordinates._globalMapping.push_back(G2L(i + _index * coordinates._points.size(), i));
	}
}

void Sphere::elements(std::vector<size_t> &bodies, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges)
{
	_block->elements(elements, 0);
	bodies = { 0, elements.size() };
}

bool Sphere::partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners)
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

void Sphere::neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours, const std::vector<Element*> &faces, const std::vector<Element*> &edges)
{
	std::vector<int> map(27, -1);
	std::vector<eslocal> cross(9 * _settings.clusters * _settings.clusters, -1);

	auto setCross = [&] (size_t offset_r, size_t offset_c, size_t row, size_t col, size_t cluster) {
		cross[(row + offset_r * _settings.clusters) * 3 * _settings.clusters + col + offset_c * _settings.clusters] = cluster;
	};

	for (size_t i = 0, index = 0, c = _settings.clusters - 1; i < _settings.clusters; i++) {
		for (size_t j = 0; j < _settings.clusters; j++, index++) {
			switch (_side) {
			case SIDE::UP:
				setCross(0, 1,     j,     i, index + static_cast<int>(SIDE::BACK)  * _settings.clusters * _settings.clusters);
				setCross(1, 0,     i, c - j, index + static_cast<int>(SIDE::LEFT)  * _settings.clusters * _settings.clusters);
				setCross(1, 1,     j,     i, index + static_cast<int>(SIDE::UP)    * _settings.clusters * _settings.clusters);
				setCross(1, 2, c - i,     j, index + static_cast<int>(SIDE::RIGHT) * _settings.clusters * _settings.clusters);
				setCross(2, 1,     j,     i, index + static_cast<int>(SIDE::FRONT) * _settings.clusters * _settings.clusters);
				break;
			case SIDE::FRONT:
				setCross(0, 1,     j,     i, index + static_cast<int>(SIDE::UP)    * _settings.clusters * _settings.clusters);
				setCross(1, 0,     j,     i, index + static_cast<int>(SIDE::LEFT)  * _settings.clusters * _settings.clusters);
				setCross(1, 1,     j,     i, index + static_cast<int>(SIDE::FRONT) * _settings.clusters * _settings.clusters);
				setCross(1, 2,     j,     i, index + static_cast<int>(SIDE::RIGHT) * _settings.clusters * _settings.clusters);
				setCross(2, 1,     j,     i, index + static_cast<int>(SIDE::DOWN)  * _settings.clusters * _settings.clusters);
				break;
			case SIDE::DOWN:
				setCross(0, 1,     j,     i, index + static_cast<int>(SIDE::FRONT) * _settings.clusters * _settings.clusters);
				setCross(1, 0, c - i,     j, index + static_cast<int>(SIDE::LEFT)  * _settings.clusters * _settings.clusters);
				setCross(1, 1,     j,     i, index + static_cast<int>(SIDE::DOWN)  * _settings.clusters * _settings.clusters);
				setCross(1, 2,     i, c - j, index + static_cast<int>(SIDE::RIGHT) * _settings.clusters * _settings.clusters);
				setCross(2, 1,     j,     i, index + static_cast<int>(SIDE::BACK)  * _settings.clusters * _settings.clusters);
				break;
			case SIDE::BACK:
				setCross(0, 1,     j,     i, index + static_cast<int>(SIDE::DOWN)  * _settings.clusters * _settings.clusters);
				setCross(1, 0, c - j, c - i, index + static_cast<int>(SIDE::LEFT)  * _settings.clusters * _settings.clusters);
				setCross(1, 1,     j,     i, index + static_cast<int>(SIDE::BACK)  * _settings.clusters * _settings.clusters);
				setCross(1, 2, c - j, c - i, index + static_cast<int>(SIDE::RIGHT) * _settings.clusters * _settings.clusters);
				setCross(2, 1,     j,     i, index + static_cast<int>(SIDE::UP)    * _settings.clusters * _settings.clusters);
				break;
			case SIDE::LEFT:
				setCross(0, 1, c - i,     j, index + static_cast<int>(SIDE::UP)    * _settings.clusters * _settings.clusters);
				setCross(1, 0, c - j, c - i, index + static_cast<int>(SIDE::BACK)  * _settings.clusters * _settings.clusters);
				setCross(1, 1,     j,     i, index + static_cast<int>(SIDE::LEFT)  * _settings.clusters * _settings.clusters);
				setCross(1, 2,     j,     i, index + static_cast<int>(SIDE::FRONT) * _settings.clusters * _settings.clusters);
				setCross(2, 1,     i, c - j, index + static_cast<int>(SIDE::DOWN)  * _settings.clusters * _settings.clusters);
				break;
			case SIDE::RIGHT:
				setCross(0, 1,     i, c - j, index + static_cast<int>(SIDE::UP)    * _settings.clusters * _settings.clusters);
				setCross(1, 0,     j,     i, index + static_cast<int>(SIDE::FRONT) * _settings.clusters * _settings.clusters);
				setCross(1, 1,     j,     i, index + static_cast<int>(SIDE::RIGHT) * _settings.clusters * _settings.clusters);
				setCross(1, 2, c - j, c - i, index + static_cast<int>(SIDE::BACK)  * _settings.clusters * _settings.clusters);
				setCross(2, 1, c - i,     j, index + static_cast<int>(SIDE::DOWN)  * _settings.clusters * _settings.clusters);
				break;
			}

		}
	}

	size_t mid = (_row + _settings.clusters) * 3 * _settings.clusters + _col + _settings.clusters;

	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			if (cross[mid + j * 3 * _settings.clusters + i] != -1) {
				map[13 + i * 3 + j] = cross[mid + j * 3 * _settings.clusters + i] + _layer * 6 * _settings.clusters * _settings.clusters;
			}
		}
	}

	if (_layer) {
		for (int i = 0; i < 9; i++) {
			if (map[i + 9] != -1) {
				map[i] = map[i + 9] - 6 * _settings.clusters * _settings.clusters;
			}
		}
	}

	if (_layer + 1 < _settings.layers) {
		for (int i = 0; i < 9; i++) {
			if (map[i + 9] != -1) {
				map[i + 18] = map[i + 9] + 6 * _settings.clusters * _settings.clusters;
			}
		}
	}

	_block->boundaries(nodes, map);
	for (size_t i = 0; i < map.size(); i++) {
		if (i != 13 && map[i] != -1) {
			neighbours.push_back(map[i]);
		}
	}
	std::sort(neighbours.begin(), neighbours.end());
	mesh.synchronizeGlobalIndices();
}

void Sphere::regions(
		std::vector<Evaluator*> &evaluators,
		std::vector<Region*> &regions,
		std::vector<Element*> &elements,
		std::vector<Element*> &faces,
		std::vector<Element*> &edges,
		std::vector<Element*> &nodes)
{
	for (auto it = _sphere.nodes.begin(); it != _sphere.nodes.end(); ++it) {
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			regions.push_back(new Region(ElementType::NODES, nodes));
		} else {
			regions.push_back(new Region(ElementType::NODES));
			BlockBorder border(it->second);
			_block->region(nodes, regions.back(), border, 0);
		}
		regions.back()->name = it->first;
	}
	for (auto it = _sphere.edges.begin(); it != _sphere.edges.end(); ++it) {
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			ESINFO(GLOBAL_ERROR) << "Implement region of all edges.";
		} else {
			regions.push_back(new Region(ElementType::EDGES));
			BlockBorder border(it->second);
			_block->region(nodes, regions.back(), border, 1);
			edges.insert(edges.end(), regions.back()->elements().begin(), regions.back()->elements().end());
		}
		regions.back()->name = it->first;
	}
	for (auto it = _sphere.faces.begin(); it != _sphere.faces.end(); ++it) {
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			ESINFO(GLOBAL_ERROR) << "Implement region of all faces.";
		} else {
			regions.push_back(new Region(ElementType::FACES));
			BlockBorder border(it->second);
			_block->region(nodes, regions.back(), border, 2);
			faces.insert(faces.end(), regions.back()->elements().begin(), regions.back()->elements().end());
		}
		regions.back()->name = it->first;
	}

	for (auto it = _sphere.elements.begin(); it != _sphere.elements.end(); ++it) {
		if (StringCompare::caseInsensitiveEq("all", it->second)) {
			regions.push_back(new Region(ElementType::ELEMENTS, elements));
		} else {
			regions.push_back(new Region(ElementType::ELEMENTS));
			BlockBorder border(it->second);
			_block->region(elements, regions.back(), border, 3);
		}
		regions.back()->name = it->first;
	}
}




