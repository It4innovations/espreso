
#include "spheregenerator.h"
#include "input/parsers/meshgenerator/meshgenerator.h"

#include "input/parsers/meshgenerator/selection/blockborder.h"
#include "input/meshbuilder.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"
#include "config/ecf/input/sphere.h"

using namespace espreso;

void SphereGenerator::generate(const SphereGeneratorConfiguration &configuration, MeshBuilder &mesh)
{
	SphereGenerator sphere(configuration);

	sphere.nodes(mesh);
	sphere.elements(mesh);
	sphere.neighbors(mesh);
	sphere.regions(configuration, mesh);
}


SphereGenerator::SphereGenerator(const SphereGeneratorConfiguration &configuration)
: GridGenerator(configuration)
{
	if (info::mpi::size % 6 != 0) {
		eslog::globalerror("Number of MPI process should be 6 x clusters x clusters for sphere generator.\n");
	}

	// clusters.x == clusters.y

	_row = ((info::mpi::rank % (_settings.clusters.x * _settings.clusters.y)) % _settings.clusters.x);
	_col = ((info::mpi::rank % (_settings.clusters.x * _settings.clusters.y)) / _settings.clusters.x);
	_layer = info::mpi::rank / (6 * _settings.clusters.x * _settings.clusters.y);

	_settings.start = Triple<esint>( _row      / _settings.clusters.x / MeshGenerator::precision,  _col      / _settings.clusters.y / MeshGenerator::precision,  _layer      / _settings.clusters.z / MeshGenerator::precision);
	_settings.end   = Triple<esint>((_row + 1) / _settings.clusters.x / MeshGenerator::precision, (_col + 1) / _settings.clusters.y / MeshGenerator::precision, (_layer + 1) / _settings.clusters.z / MeshGenerator::precision);

	switch ((info::mpi::rank / (_settings.clusters.x * _settings.clusters.y)) % 6) {
	case 0:
		_side = SIDE::UP;
		_settings.projection = Triple<Expression>(
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(x * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(y * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4))"                           , { "x", "y", "z" }));
		break;
	case 1:
		_side = SIDE::FRONT;
		_settings.projection = Triple<Expression>(
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4))"                           , { "x", "y", "z" }),
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(y * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(x * pi / 2 - pi / 4)", { "x", "y", "z" }));;
		break;
	case 2:
		_side = SIDE::DOWN;
		_settings.projection = Triple<Expression>(
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(x * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(y * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4))"                           , { "x", "y", "z" }));
		break;
	case 3:
		_side = SIDE::BACK;
		_settings.projection = Triple<Expression>(
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4))"                           , { "x", "y", "z" }),
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(y * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(x * pi / 2 - pi / 4)", { "x", "y", "z" }));;
		break;
	case 4:
		_side = SIDE::LEFT;
		_settings.projection = Triple<Expression>(
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(y * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4))"                           , { "x", "y", "z" }),
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(x * pi / 2 - pi / 4)", { "x", "y", "z" }));
		break;
	case 5:
		_side = SIDE::RIGHT;
		_settings.projection = Triple<Expression>(
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(y * pi / 2 - pi / 4)", { "x", "y", "z" }),
			Expression(" (5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4))"                           , { "x", "y", "z" }),
			Expression("-(5 + 5 * z) / sqrt(1 + tan(x * pi / 2 - pi / 4) * tan(x * pi / 2 - pi / 4) + tan(y * pi / 2 - pi / 4) * tan(y * pi / 2 - pi / 4)) * tan(x * pi / 2 - pi / 4)", { "x", "y", "z" }));
		break;
	}
}

void SphereGenerator::nodes(MeshBuilder &mesh)
{
	GridGenerator::nodes(mesh);
	for (auto n = mesh.nIDs.begin(); n != mesh.nIDs.end(); ++n) {
		*n += info::mpi::rank * mesh.nIDs.size();
	}
}

void SphereGenerator::neighbors(MeshBuilder &mesh)
{
	std::vector<int> map(27, -1);
	std::vector<esint> cross(9 * _settings.clusters.x * _settings.clusters.y, -1);

	auto setCross = [&] (size_t offset_r, size_t offset_c, size_t row, size_t col, size_t cluster) {
		cross[(row + offset_r * _settings.clusters.x) * 3 * _settings.clusters.y + col + offset_c * _settings.clusters.x] = cluster;
	};

	for (size_t i = 0, index = 0, c = _settings.clusters.x - 1; i < _settings.clusters.x; i++) {
		for (size_t j = 0; j < _settings.clusters.y; j++, index++) {
			switch (_side) {
			case SIDE::UP:
				setCross(0, 1,     j,     i, index + static_cast<int>(SIDE::BACK)  * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 0,     i, c - j, index + static_cast<int>(SIDE::LEFT)  * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 1,     j,     i, index + static_cast<int>(SIDE::UP)    * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 2, c - i,     j, index + static_cast<int>(SIDE::RIGHT) * _settings.clusters.x * _settings.clusters.y);
				setCross(2, 1,     j,     i, index + static_cast<int>(SIDE::FRONT) * _settings.clusters.x * _settings.clusters.y);
				break;
			case SIDE::FRONT:
				setCross(0, 1,     j,     i, index + static_cast<int>(SIDE::UP)    * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 0,     j,     i, index + static_cast<int>(SIDE::LEFT)  * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 1,     j,     i, index + static_cast<int>(SIDE::FRONT) * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 2,     j,     i, index + static_cast<int>(SIDE::RIGHT) * _settings.clusters.x * _settings.clusters.y);
				setCross(2, 1,     j,     i, index + static_cast<int>(SIDE::DOWN)  * _settings.clusters.x * _settings.clusters.y);
				break;
			case SIDE::DOWN:
				setCross(0, 1,     j,     i, index + static_cast<int>(SIDE::FRONT) * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 0, c - i,     j, index + static_cast<int>(SIDE::LEFT)  * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 1,     j,     i, index + static_cast<int>(SIDE::DOWN)  * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 2,     i, c - j, index + static_cast<int>(SIDE::RIGHT) * _settings.clusters.x * _settings.clusters.y);
				setCross(2, 1,     j,     i, index + static_cast<int>(SIDE::BACK)  * _settings.clusters.x * _settings.clusters.y);
				break;
			case SIDE::BACK:
				setCross(0, 1,     j,     i, index + static_cast<int>(SIDE::DOWN)  * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 0, c - j, c - i, index + static_cast<int>(SIDE::LEFT)  * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 1,     j,     i, index + static_cast<int>(SIDE::BACK)  * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 2, c - j, c - i, index + static_cast<int>(SIDE::RIGHT) * _settings.clusters.x * _settings.clusters.y);
				setCross(2, 1,     j,     i, index + static_cast<int>(SIDE::UP)    * _settings.clusters.x * _settings.clusters.y);
				break;
			case SIDE::LEFT:
				setCross(0, 1, c - i,     j, index + static_cast<int>(SIDE::UP)    * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 0, c - j, c - i, index + static_cast<int>(SIDE::BACK)  * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 1,     j,     i, index + static_cast<int>(SIDE::LEFT)  * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 2,     j,     i, index + static_cast<int>(SIDE::FRONT) * _settings.clusters.x * _settings.clusters.y);
				setCross(2, 1,     i, c - j, index + static_cast<int>(SIDE::DOWN)  * _settings.clusters.x * _settings.clusters.y);
				break;
			case SIDE::RIGHT:
				setCross(0, 1,     i, c - j, index + static_cast<int>(SIDE::UP)    * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 0,     j,     i, index + static_cast<int>(SIDE::FRONT) * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 1,     j,     i, index + static_cast<int>(SIDE::RIGHT) * _settings.clusters.x * _settings.clusters.y);
				setCross(1, 2, c - j, c - i, index + static_cast<int>(SIDE::BACK)  * _settings.clusters.x * _settings.clusters.y);
				setCross(2, 1, c - i,     j, index + static_cast<int>(SIDE::DOWN)  * _settings.clusters.x * _settings.clusters.y);
				break;
			}

		}
	}

	size_t mid = (_row + _settings.clusters.x) * 3 * _settings.clusters.y + _col + _settings.clusters.x;

	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			if (cross[mid + j * 3 * _settings.clusters.x + i] != -1) {
				if (cross[mid + j * 3 * _settings.clusters.x + i] + _layer * 6 * _settings.clusters.x * _settings.clusters.y < (size_t)info::mpi::size) {
					map[13 + i * 3 + j] = cross[mid + j * 3 * _settings.clusters.x + i] + _layer * 6 * _settings.clusters.x * _settings.clusters.y;
				}
			}
		}
	}

	if (_layer) {
		for (int i = 0; i < 9; i++) {
			if (map[i + 9] != -1) {
				map[i] = map[i + 9] - 6 * _settings.clusters.x * _settings.clusters.y;
			}
		}
	}

	if (_layer + 1 < _settings.clusters.z) {
		for (int i = 0; i < 9; i++) {
			if (map[i + 9] != -1) {
				map[i + 18] = map[i + 9] + 6 * _settings.clusters.x * _settings.clusters.y;
			}
		}
	}

	_block.neighbors(map, mesh);
}

void SphereGenerator::regions(const SphereGeneratorConfiguration &configuration, MeshBuilder &mesh)
{
	for (auto it = configuration.nodes.begin(); it != configuration.nodes.end(); ++it) {
		_block.nodesRegion(BlockBorder(it->second, configuration), mesh.nregions[it->first]);
	}
	for (auto it = configuration.edges.begin(); it != configuration.edges.end(); ++it) {
		_block.edgesRegion(BlockBorder(it->second, configuration), mesh, mesh.eregions[it->first]);
	}
	for (auto it = configuration.faces.begin(); it != configuration.faces.end(); ++it) {
		_block.facesRegion(BlockBorder(it->second, configuration), mesh, mesh.eregions[it->first]);
	}
	for (auto it = configuration.elements.begin(); it != configuration.elements.end(); ++it) {
		_block.elementsRegion(BlockBorder(it->second, configuration), mesh.eregions[it->first]);
	}
}


