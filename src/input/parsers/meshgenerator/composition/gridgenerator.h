
#ifndef SRC_INPUT_GENERATOR_COMPOSITION_GRIDGENERATOR_H_
#define SRC_INPUT_GENERATOR_COMPOSITION_GRIDGENERATOR_H_

#include "input/parsers/meshgenerator/primitives/gridsettings.h"

namespace espreso {

struct MeshBuilder;
struct GridGeneratorConfiguration;
struct SphereGeneratorConfiguration;

class GridGenerator {

    friend class GridTowerGenerator;
    friend class SphereGenerator;

public:
    static void generate(const GridGeneratorConfiguration &configuration, MeshBuilder &mesh);

protected:
    GridGenerator(const GridGeneratorConfiguration &configuration);
    GridGenerator(const SphereGeneratorConfiguration &configuration);
    virtual ~GridGenerator() {}

    void init();
    virtual void nodes(MeshBuilder &mesh);
    void elements(MeshBuilder &mesh);
    void neighbors(MeshBuilder &mesh);
    void regions(const GridGeneratorConfiguration &configuration, MeshBuilder &mesh);

    GridSettings _settings;
    BlockGenerator _block;
    Triple<size_t> _clusterOffset;
    std::vector<int> _clusterIndices;
};

}


#endif /* SRC_INPUT_GENERATOR_COMPOSITION_GRIDGENERATOR_H_ */
