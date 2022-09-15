
#include "gridgenerator.h"
#include "input/parsers/meshgenerator/meshgenerator.h"

#include "input/parsers/meshgenerator/selection/blockborder.h"
#include "input/meshbuilder.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"
#include "config/ecf/input/grid.h"
#include "config/ecf/input/sphere.h"
#include "basis/utilities/parser.h"

using namespace espreso;

void GridGenerator::generate(const GridGeneratorConfiguration &configuration, MeshBuilder &mesh)
{
	eslog::startln("GRID GENERATOR: STARTED", "GRID GENERATOR");
	GridGenerator grid(configuration);

	grid.init();
	eslog::checkpointln("GRID GENERATOR: CONFIGURED");
	grid.nodes(mesh);
	eslog::checkpointln("GRID GENERATOR: NODES GENERATED");
	grid.elements(mesh);
	eslog::checkpointln("GRID GENERATOR: ELEMENTS GENERATED");
	grid.neighbors(mesh);
	eslog::checkpointln("GRID GENERATOR: NEIGHBORS COMPUTED");
	grid.regions(configuration, mesh);
	eslog::endln("GRID GENERATOR: REGIONS CREATED");
}

GridGenerator::GridGenerator(const GridGeneratorConfiguration &configuration)
: _settings(configuration), _block(configuration, _settings)
{

}

GridGenerator::GridGenerator(const SphereGeneratorConfiguration &configuration)
: _settings(configuration), _block(configuration, _settings)
{

}

void GridGenerator::init()
{
	Triple<size_t> clusters = _settings.blocks * _settings.clusters;

	int cluster = 0;
	Triple<size_t> offset;
	for (offset.z = 0; offset.z < clusters.z; offset.z++) {
		for (offset.y = 0; offset.y < clusters.y; offset.y++) {
			for (offset.x = 0; offset.x < clusters.x; offset.x++) {

				if (_settings.nonempty[cluster]) {
					_clusterIndices.push_back(cluster);
				} else {
					_clusterIndices.push_back(-1);
					continue;
				}

				if (cluster++ == info::mpi::rank) {
					_clusterOffset = offset;
					Triple<int> start = _settings.start;
					_settings.start = start + ((_settings.end - start) / (Triple<double>)_settings.clusters * offset).round();
					_settings.end   = start + ((_settings.end - start) / (Triple<double>)_settings.clusters * (offset + 1)).round();
				}
			}
		}
	}
	if (cluster != info::mpi::size) {
		eslog::globalerror("Incorrect number of MPI processes (%d). Should be %d\n", info::mpi::size, cluster);
	}
}

void GridGenerator::nodes(MeshBuilder &mesh)
{
	_block.coordinates(mesh);
	mesh.nIDs.reserve(mesh.coordinates.size());

	Triple<size_t> cnodes = _settings.domains * _settings.elements * (Triple<size_t>(_block.element()->subnodes) - 1);
	Triple<size_t> gnodes = _settings.blocks * _settings.clusters * cnodes;
	Triple<size_t> coffset = _clusterOffset * cnodes;
	Triple<size_t> size = (gnodes + 1).toSize();

	Triple<size_t> offset;
	for (offset.z = 0; offset.z <= cnodes.z; offset.z++) {
		for (offset.y = 0; offset.y <= cnodes.y; offset.y++) {
			for (offset.x = 0; offset.x <= cnodes.x; offset.x++) {
				mesh.nIDs.push_back(((coffset + offset) * size).sum());
			}
		}
	}
}

void GridGenerator::elements(MeshBuilder &mesh)
{
	_block.elements(mesh);

	mesh.material.resize(mesh.etype.size(), 0);
	mesh.body.resize(mesh.etype.size(), _settings.body);
}

void GridGenerator::neighbors(MeshBuilder &mesh)
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
					map[index] = _clusterIndices[((coffset + offset) * size).sum()];
				}
			}
		}
	}

	for (size_t i = 0; i < map.size(); i++) {
		if (map[i] >= info::mpi::size) {
			map[i] = -1;
		}
	}

	_block.neighbors(map, mesh);
}

void GridGenerator::regions(const GridGeneratorConfiguration &configuration, MeshBuilder &mesh)
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
		if (StringCompare::caseInsensitiveEq("chessboard_white", it->second)) {
			_block.pattern(_clusterOffset, _settings.blocks * _settings.clusters, mesh.eregions[it->first], Pattern::CHESSBOARD_WHITE, configuration.chessboard_size);
		} else if (StringCompare::caseInsensitiveEq("chessboard_black", it->second)) {
			_block.pattern(_clusterOffset, _settings.blocks * _settings.clusters, mesh.eregions[it->first], Pattern::CHESSBOARD_BLACK, configuration.chessboard_size);
		} else {
			_block.elementsRegion(BlockBorder(it->second, configuration), mesh.eregions[it->first]);
		}
	}
}
