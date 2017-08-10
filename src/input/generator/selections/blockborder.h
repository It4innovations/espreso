
#ifndef SRC_INPUT_GENERATOR_SELECTIONS_BLOCKBORDER_H_
#define SRC_INPUT_GENERATOR_SELECTIONS_BLOCKBORDER_H_

#include "../primitives/triple.h"
#include "../elements/element.h"

namespace espreso {
namespace input {

struct BlockSetting;

struct BlockBorder {
	BlockBorder(const std::string &interval);

	Triple<esglobal> start, end;
	Triple<bool> excludeStart, excludeEnd;

	size_t dimension() const;
	BlockBorder intersect(const BlockSetting &block) const;

	CubeFace getFace(const BlockSetting &block) const;
	CubeEdge getEdge(const BlockSetting &block) const;
};

}
}


#endif /* SRC_INPUT_GENERATOR_SELECTIONS_BLOCKBORDER_H_ */

