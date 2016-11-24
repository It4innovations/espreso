
#ifndef SRC_INPUT_GENERATOR_PRIMITIVES_BLOCK_H_
#define SRC_INPUT_GENERATOR_PRIMITIVES_BLOCK_H_

#include "clustergenerator.h"
#include "triple.h"
#include "../selections/blockborder.h"

namespace espreso {
namespace input {

struct BlockSetting {
	Triple<size_t> domains;
	Triple<size_t> elements;

	Triple<double> start, end;
};

class BlockGenerator: public ClusterGenerator {

public:
	virtual void region(const std::vector<Element*> &nodes, Region &region, const BlockBorder &border, size_t dimension) =0;
	virtual ~BlockGenerator() {};

protected:
	BlockGenerator(Mesh &mesh, const BlockSetting &block): ClusterGenerator(mesh), block(block) {};

	BlockSetting block;
};

template <class TElement>
class Block: public BlockGenerator {

public:
	Block(Mesh &mesh, const BlockSetting &block): BlockGenerator(mesh, block) {};

	void points(std::vector<Point> &points);
	void elements(std::vector<Element*> &elements);
	void boundaries(std::vector<Element*> &nodes, const std::vector<int> &neighbours);
	void uniformPartition(std::vector<eslocal> &partsPtrs, size_t subdomains);

	void region(const std::vector<Element*> &nodes, Region &region, const BlockBorder &border, size_t dimension);

private:
	void forEachElement(const Triple<size_t> &start, const Triple<size_t> &end, std::function<void(std::vector<eslocal> &indices)> operation);
	void forEachElement(const Triple<size_t> &start, const Triple<size_t> &end, std::function<void(std::vector<eslocal> &indices)> operation, std::function<void(Triple<size_t> &offset)> restriction);
};

}
}

#include "block.hpp"



#endif /* SRC_INPUT_GENERATOR_PRIMITIVES_BLOCK_H_ */
