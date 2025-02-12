
#include "gridsettings.h"

#include "config/ecf/input/grid.h"
#include "config/ecf/input/sphere.h"
#include "esinfo/eslog.h"

using namespace espreso;

GridSettings::GridSettings(const GridGeneratorConfiguration &configuration)
: BlockSettings(configuration)
{
    blocks   = Triple<size_t>(configuration.blocks_x, configuration.blocks_y, configuration.blocks_z);
    clusters = Triple<size_t>(configuration.clusters_x, configuration.clusters_y, configuration.clusters_z);
    nonempty.resize((blocks * clusters).mul(), true);

    for (auto it = configuration.blocks.begin(); it != configuration.blocks.end(); ++it) {
        if (it->first >= nonempty.size()) {
            eslog::globalerror("Block index is out of range.\n");
        }
        nonempty[it->first] = it->second;
    }
    body = 0;
}

GridSettings::GridSettings(const SphereGeneratorConfiguration &configuration)
: BlockSettings(configuration)
{
    blocks   = Triple<size_t>(1, 1, 1);
    clusters = Triple<size_t>(configuration.clusters, configuration.clusters, configuration.layers);
    nonempty.resize((blocks * clusters).mul(), true);
    body = 0;
}


