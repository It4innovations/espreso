
#include "gridsetgenerator.h"
#include "input/parsers/meshgenerator/meshgenerator.h"

#include "input/parsers/meshgenerator/elements/element.h"
#include "input/meshbuilder.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"
#include "config/ecf/input/gridset.h"
#include "wrappers/mpi/communication.h"

using namespace espreso;

void GridSetGenerator::generate(const GridSetGeneratorConfiguration &configuration, MeshBuilder &mesh)
{
	eslog::startln("GRID SET GENERATOR: STARTED", "GRID SET GENERATOR");
	GridSetGenerator set(configuration);

	set.init(configuration);
	eslog::checkpointln("GRID SET GENERATOR: CONFIGURED");
	set.nodes(mesh);
	eslog::checkpointln("GRID SET GENERATOR: NODES GENERATED");
	set.elements(mesh);
	eslog::checkpointln("GRID SET GENERATOR: ELEMENTS GENERATED");
	set.neighbors(mesh);
	eslog::checkpointln("GRID SET GENERATOR: NEIGHBORS COMPUTED");
	set.regions(configuration, mesh);
	eslog::endln("GRID SET GENERATOR: REGIONS CREATED");
}

esint GridSetGenerator::grids(const GridSetGeneratorConfiguration &configuration)
{
	int clusters = 0;
	for (auto grid = configuration.grids.begin(); grid != configuration.grids.end(); ++grid) {
		size_t nclusters = grid->second.clusters_x * grid->second.clusters_y * grid->second.clusters_z;
		clusters += nclusters;
		for (auto it = grid->second.blocks.begin(); it != grid->second.blocks.end(); ++it) {
			if (it->first >= nclusters) {
				eslog::globalerror("Block index is out of range.\n");
			}
			if (!it->second) {
				--clusters;
			}
		}
	}

	return clusters;
}

esint GridSetGenerator::gridIndex(const GridSetGeneratorConfiguration &configuration)
{
	int clusters = 0;
	for (auto grid = configuration.grids.begin(); grid != configuration.grids.end(); ++grid) {
		size_t maxclusters = grid->second.clusters_x * grid->second.clusters_y * grid->second.clusters_z;
		int gclusters = maxclusters;
		for (auto it = grid->second.blocks.begin(); it != grid->second.blocks.end(); ++it) {
			if (it->first >= maxclusters) {
				eslog::globalerror("Block index is out of range.\n");
			}
			if (!it->second) {
				--gclusters;
			}
		}

		if (clusters <= info::mpi::rank && info::mpi::rank < clusters + gclusters) {
			return grid->first;
		}
		clusters += gclusters;
	}

	return 0;
}


GridSetGenerator::GridSetGenerator(const GridSetGeneratorConfiguration &configuration)
: GridGenerator(configuration.grids.at(gridIndex(configuration))), _gridNodeOffset(0)
{

}

void GridSetGenerator::init(const GridSetGeneratorConfiguration &configuration)
{
	_gridIndex = gridIndex(configuration);
	Triple<size_t> clusters = _settings.blocks * _settings.clusters;

	int cluster = 0;
	size_t nodes = 0, noffset = 0;

	Triple<esint> lenght;
	for (auto grid = configuration.grids.begin(); grid != configuration.grids.end(); ++grid) {
		if (grid->first == _gridIndex) {
			Triple<size_t> offset;
			int begin = cluster;
			for (offset.z = 0; offset.z < clusters.z; offset.z++) {
				for (offset.y = 0; offset.y < clusters.y; offset.y++) {
					for (offset.x = 0; offset.x < clusters.x; offset.x++) {
						if (_settings.nonempty[cluster - begin]) {
							_clusterIndices.push_back(cluster);
						} else {
							_clusterIndices.push_back(-1);
							continue;
						}

						if (cluster++ == info::mpi::rank) {
							_clusterOffset = offset;
							_settings.start = _settings.start + ((_settings.end - _settings.start) / (Triple<double>)_settings.clusters * offset).round();
							_settings.end   = _settings.start + ((_settings.end - _settings.start) / (Triple<double>)_settings.clusters * (offset + 1)).round();
						}
					}
				}
			}
			nodes = (_settings.blocks * _settings.clusters * _settings.domains * _settings.elements * (Triple<size_t>(_block.element()->subnodes) - 1) + 1).mul();
		}
		int maxcluster = cluster;
		Communication::allReduce(&maxcluster, &cluster, 1, MPI_INT, MPI_MAX);
		noffset += nodes;
		Communication::allReduce(&noffset, &nodes, 1, MPITools::getType<size_t>().mpitype, MPI_MAX);
		noffset = nodes;
		if (grid->first + 1 == _gridIndex) {
			_gridNodeOffset = noffset;
		}
		lenght += _settings.end - _settings.start;
	}

	_settings.body = _gridIndex;

	if (grids(configuration) != info::mpi::size) {
		eslog::globalerror("Incorrect number of MPI processes (%d). Should be %d\n", info::mpi::size, grids(configuration));
	}
}

void GridSetGenerator::nodes(MeshBuilder &mesh)
{
	GridGenerator::nodes(mesh);
	for (auto n = mesh.nIDs.begin(); n != mesh.nIDs.end(); ++n) {
		*n += _gridNodeOffset;
	}
}

void GridSetGenerator::regions(const GridSetGeneratorConfiguration &configuration, MeshBuilder &mesh)
{
	for (auto it = configuration.grids.begin(); it != configuration.grids.end(); ++it) {
		GridGenerator::regions(it->second, mesh);
	}
}
