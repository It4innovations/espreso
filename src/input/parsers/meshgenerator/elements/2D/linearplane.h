
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_2D_LINEARPLANE_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_2D_LINEARPLANE_H_

#include "input/parsers/meshgenerator/elements/element.h"

namespace espreso {

struct LinearPlaneGenerator: public ElementGenerator {

	LinearPlaneGenerator();

	void pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeEdge edge) const;
	void pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const;

	void pushEdge(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeEdge edge) const;
	void pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const;
};

}


#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_2D_LINEARPLANE_H_ */
