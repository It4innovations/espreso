
#ifndef SRC_CONFIG_ECF_INPUT_GRID_H_
#define SRC_CONFIG_ECF_INPUT_GRID_H_

#include "block.h"

#include <map>

namespace espreso {

struct GridGeneratorConfiguration: public BlockGeneratorConfiguration {

    size_t blocks_x, blocks_y, blocks_z;
    size_t clusters_x, clusters_y, clusters_z;

    std::map<size_t, bool> blocks;
    std::map<size_t, size_t> noncontinuous;

    std::map<std::string, std::string> nodes, edges, faces, elements;

    size_t chessboard_size;
    size_t nonuniform_nparts;

    GridGeneratorConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_GRID_H_ */
