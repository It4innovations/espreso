
#include "gridtowergenerator.h"
#include "input/parsers/meshgenerator/meshgenerator.h"

#include "input/parsers/meshgenerator/elements/element.h"
#include "input/meshbuilder.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"
#include "config/ecf/input/gridtower.h"
#include "wrappers/mpi/communication.h"

using namespace espreso;

void GridTowerGenerator::generate(const GridTowerGeneratorConfiguration &configuration, MeshBuilder &mesh)
{
    eslog::startln("GRID TOWER GENERATOR: STARTED", "GRID TOWER GENERATOR");
    GridTowerGenerator tower(configuration);

    tower.init(configuration);
    eslog::checkpointln("GRID TOWER GENERATOR: CONFIGURED");
    tower.nodes(mesh);
    eslog::checkpointln("GRID TOWER GENERATOR: NODES GENERATED");
    tower.elements(mesh);
    eslog::checkpointln("GRID TOWER GENERATOR: ELEMENTS GENERATED");
    tower.neighbors(configuration, mesh);
    eslog::checkpointln("GRID TOWER GENERATOR: NEIGHBORS COMPUTED");
    tower.regions(configuration, mesh);
    eslog::endln("GRID TOWER GENERATOR: REGIONS CREATED");
}

esint GridTowerGenerator::grids(const GridTowerGeneratorConfiguration &configuration)
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

esint GridTowerGenerator::gridIndex(const GridTowerGeneratorConfiguration &configuration)
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


GridTowerGenerator::GridTowerGenerator(const GridTowerGeneratorConfiguration &configuration)
: GridGenerator(configuration.grids.at(gridIndex(configuration))), _gridNodeOffset(0)
{

}

void GridTowerGenerator::init(const GridTowerGeneratorConfiguration &configuration)
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
                            Triple<long> start = _settings.start;
                            _settings.start = start + ((_settings.end - start) / (Triple<double>)_settings.clusters * offset).round();
                            _settings.end   = start + ((_settings.end - start) / (Triple<double>)_settings.clusters * (offset + 1)).round();
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
        if (configuration.composition == GridTowerGeneratorConfiguration::COMPOSITION::GLUED) {
            if (grid->first == _gridIndex) {
                switch (configuration.direction) {
                case GridTowerGeneratorConfiguration::DIRECTION::X:
                        _settings.end.x += lenght.x;
                        _settings.start.x += lenght.x;
                    break;
                case GridTowerGeneratorConfiguration::DIRECTION::Y:
                        _settings.end.y += lenght.y;
                        _settings.start.y += lenght.y;;
                    break;
                case GridTowerGeneratorConfiguration::DIRECTION::Z:
                        _settings.end.z += lenght.z;
                        _settings.start.z += lenght.z;
                    break;
                }
            }
            lenght += _settings.end - _settings.start;
        }
    }

    _settings.body = _gridIndex;

    if (grids(configuration) != info::mpi::size) {
        eslog::globalerror("Incorrect number of MPI processes (%d). Should be %d\n", info::mpi::size, grids(configuration));
    }
}

void GridTowerGenerator::nodes(MeshBuilder &mesh)
{
    GridGenerator::nodes(mesh);
    for (auto n = mesh.nIDs.begin(); n != mesh.nIDs.end(); ++n) {
        *n += _gridNodeOffset;
    }
}

void GridTowerGenerator::neighbors(const GridTowerGeneratorConfiguration &configuration, MeshBuilder &mesh)
{
    auto middle = configuration.grids.find(_gridIndex);
    auto lower  = configuration.grids.end();
    auto upper  = configuration.grids.end();

    auto lowerUpper = [&] (size_t offset, size_t max) {
        if (offset == 0) {
            lower  = configuration.grids.find(_gridIndex - 1);
        }
        if (offset == max) {
            upper  = configuration.grids.find(_gridIndex + 1);
        }
    };

    switch (configuration.direction) {
    case GridTowerGeneratorConfiguration::DIRECTION::X:
        lowerUpper(_clusterOffset.x, _settings.clusters.x - 1);
        break;
    case GridTowerGeneratorConfiguration::DIRECTION::Y:
        lowerUpper(_clusterOffset.y, _settings.clusters.y - 1);
        break;
    case GridTowerGeneratorConfiguration::DIRECTION::Z:
        lowerUpper(_clusterOffset.z, _settings.clusters.z - 1);
        break;
    }

    if (lower != configuration.grids.end()) {
        if (lower->second.blocks_z != middle->second.blocks_z || lower->second.blocks_y != middle->second.blocks_y || lower->second.blocks_x != middle->second.blocks_x) {
            lower = configuration.grids.end();
        }
        if (lower->second.clusters_z != middle->second.clusters_z || lower->second.clusters_y != middle->second.clusters_y || lower->second.clusters_x != middle->second.clusters_x) {
            lower = configuration.grids.end();
        }
        if (lower->second.domains_z != middle->second.domains_z || lower->second.domains_y != middle->second.domains_y || lower->second.domains_x != middle->second.domains_x) {
            lower = configuration.grids.end();
        }
        if (lower->second.elements_z != middle->second.elements_z || lower->second.elements_y != middle->second.elements_y || lower->second.elements_x != middle->second.elements_x) {
            lower = configuration.grids.end();
        }
    }
    if (upper != configuration.grids.end()) {
        if (upper->second.blocks_z != middle->second.blocks_z || upper->second.blocks_y != middle->second.blocks_y || upper->second.blocks_x != middle->second.blocks_x) {
            upper = configuration.grids.end();
        }
        if (upper->second.clusters_z != middle->second.clusters_z || upper->second.clusters_y != middle->second.clusters_y || upper->second.clusters_x != middle->second.clusters_x) {
            upper = configuration.grids.end();
        }
        if (upper->second.domains_z != middle->second.domains_z || upper->second.domains_y != middle->second.domains_y || upper->second.domains_x != middle->second.domains_x) {
            upper = configuration.grids.end();
        }
        if (upper->second.elements_z != middle->second.elements_z || upper->second.elements_y != middle->second.elements_y || upper->second.elements_x != middle->second.elements_x) {
            upper = configuration.grids.end();
        }
    }

    // temporary disable matching grids
    lower = configuration.grids.end();
    upper = configuration.grids.end();

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


    size_t cx = middle->second.blocks_x * middle->second.clusters_x;
    size_t cy = middle->second.blocks_y * middle->second.clusters_y;
    size_t cz = middle->second.blocks_z * middle->second.clusters_z;

    for (offset.z = -1, index = 0; offset.z <= 1; offset.z++) {
        for (offset.y = -1; offset.y <= 1; offset.y++) {
            for (offset.x = -1; offset.x <= 1; offset.x++, index++) {
                switch (configuration.direction) {
                case GridTowerGeneratorConfiguration::DIRECTION::X:
                    if (offset.x == 1 && map[index - 1] != -1 && upper != configuration.grids.end()) {
                        map[index] = map[index - 1] + cx * cy * (cz - 1) + cx * (cy - 1) + 1;
                    }
                    if (offset.x == -1 && map[index + 1] != -1 && lower != configuration.grids.end()) {
                        map[index] = map[index + 1] - cx * cy * (cz - 1) - cx * (cy - 1) - 1;
                    }
                    break;
                case GridTowerGeneratorConfiguration::DIRECTION::Y:
                    if (offset.y == 1 && map[index - 3] != -1 && upper != configuration.grids.end()) {
                        map[index] = map[index - 3] + cx * cy * (cz - 1) + cx;
                    }
                    if (offset.y == -1 && map[index + 3] != -1 && lower != configuration.grids.end()) {
                        map[index] = map[index + 3] - cx * cy * (cz - 1) - cx;
                    }
                    break;
                case GridTowerGeneratorConfiguration::DIRECTION::Z:
                    if (offset.z == 1 && map[index - 9] != -1 && upper != configuration.grids.end()) {
                        map[index] = map[index - 9] + cx * cy;
                    }
                    if (offset.z == -1 && map[index + 9] != -1 && lower != configuration.grids.end()) {
                        map[index] = map[index + 9] - cx * cy;
                    }
                    break;
                }
            }
        }
    }

    _block.neighbors(map, mesh);
}

void GridTowerGenerator::regions(const GridTowerGeneratorConfiguration &configuration, MeshBuilder &mesh)
{
    for (auto it = configuration.grids.begin(); it != configuration.grids.end(); ++it) {
        GridGenerator::regions(it->second, mesh);
    }
}
