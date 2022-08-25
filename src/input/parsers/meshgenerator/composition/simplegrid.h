
#ifndef SRC_INPUT_PARSERS_MESHGENERATOR_COMPOSITION_SIMPLEGRID_H_
#define SRC_INPUT_PARSERS_MESHGENERATOR_COMPOSITION_SIMPLEGRID_H_

#include "generator.h"
#include "input/parsers/meshgenerator/primitives/gridsettings.h"

namespace espreso {

class SimpleGridGenerator: public Generator {

public:
	SimpleGridGenerator(const GridGeneratorConfiguration &configuration);

	void nodes(NodesDomain &nodes);
	void elements(Elements &elements);
	void neighbors(Domain &domain);
	void regions();

protected:
	GridSettings _settings;
	BlockGenerator _block;
	Triple<size_t> _clusterOffset;
	std::vector<int> _clusterIndices;
};

}



#endif /* SRC_INPUT_PARSERS_MESHGENERATOR_COMPOSITION_SIMPLEGRID_H_ */
