
#include "simplegrid.h"

#include "config/ecf/input/grid.h"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"
#include "input/parsers/meshgenerator/selection/blockborder.h"

using namespace espreso;

SimpleGridGenerator::SimpleGridGenerator(const GridGeneratorConfiguration &configuration)
: _settings(configuration), _block(configuration, _settings)
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

void SimpleGridGenerator::nodes(NodesDomain &nodes)
{
	_block.coordinates(nodes);
	nodes.ids.reserve(nodes.coordinates.size());

	Triple<size_t> cnodes = _settings.domains * _settings.elements * (Triple<size_t>(_block.element()->subnodes) - 1);
	Triple<size_t> gnodes = _settings.blocks * _settings.clusters * cnodes;
	Triple<size_t> coffset = _clusterOffset * cnodes;
	Triple<size_t> size = (gnodes + 1).toSize();

	Triple<size_t> offset;
	for (offset.z = 0; offset.z <= cnodes.z; offset.z++) {
		for (offset.y = 0; offset.y <= cnodes.y; offset.y++) {
			for (offset.x = 0; offset.x <= cnodes.x; offset.x++) {
				nodes.ids.push_back(((coffset + offset) * size).sum());
			}
		}
	}
}

void SimpleGridGenerator::elements(Elements &elements)
{
	_block.elements(elements);
}

void SimpleGridGenerator::neighbors(Domain &domain)
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

	_block.neighbors(map, domain);
}

void SimpleGridGenerator::regions()
{

}
