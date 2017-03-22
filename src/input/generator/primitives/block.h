
#ifndef SRC_INPUT_GENERATOR_PRIMITIVES_BLOCK_H_
#define SRC_INPUT_GENERATOR_PRIMITIVES_BLOCK_H_

#include "clustergenerator.h"

#include "triple.h"
#include "blocksetting.h"
#include "../selections/blockborder.h"

#include "../../../mesh/elements/element.h"

namespace espreso {
namespace input {


class BlockGenerator: public ClusterGenerator {

public:
	virtual void region(const std::vector<Element*> &elements, Region *region, const BlockBorder &border, size_t dimension) =0;
	virtual ~BlockGenerator() {};

	BlockSetting block;
protected:
	BlockGenerator(Mesh &mesh, const BlockSetting &block): block(block) {};
};

template <class TElement>
class Block: public BlockGenerator {

public:
	Block(Mesh &mesh, const BlockSetting &block): BlockGenerator(mesh, block) {};

	void points(std::vector<Point> &points);
	void elements(std::vector<Element*> &elements, size_t body);
	void boundaries(std::vector<Element*> &nodes, const std::vector<int> &neighbours);
	void uniformPartition(std::vector<eslocal> &partsPtrs, size_t subdomains);
	void uniformFixPoints(const std::vector<Element*> &nodes, std::vector<std::vector<Element*> > &fixPoints);
	void uniformCorners(const std::vector<Element*> &nodes, std::vector<Element*> &corners, size_t number, bool point, bool edge, bool face);

	void region(const std::vector<Element*> &elements, Region *region, const BlockBorder &border, size_t dimension);

private:
	void forEachElement(const Triple<size_t> &start, const Triple<size_t> &end, std::function<void(std::vector<eslocal> &indices)> operation);
	void forEachElement(const Triple<size_t> &start, const Triple<size_t> &end, std::function<void(std::vector<eslocal> &indices)> operation, std::function<void(Triple<size_t> &offset)> restriction);
};

}
}

#include "block.hpp"



#endif /* SRC_INPUT_GENERATOR_PRIMITIVES_BLOCK_H_ */
