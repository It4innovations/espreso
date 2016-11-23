
#ifndef SRC_INPUT_GENERATOR_SELECTIONS_BLOCKBORDER_H_
#define SRC_INPUT_GENERATOR_SELECTIONS_BLOCKBORDER_H_

#include "../triple.h"
#include "../../meshgenerator/elements/element.h"

namespace espreso {
namespace input {

struct BlockSetting;

struct BlockBorder {
	BlockBorder(const std::string &interval);

	double epsilon = 0.001;
	Triple<double> start, end;
	Triple<bool> excludeStart, excludeEnd;

	size_t dimension() const;
	bool intersect(const BlockSetting &block) const;

	CubeFace getFace(const BlockSetting &block) const;
	CubeEdge getEdge(const BlockSetting &block) const;
};

}
}


#endif /* SRC_INPUT_GENERATOR_SELECTIONS_BLOCKBORDER_H_ */

