
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_ELEMENT_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_ELEMENT_H_

#include "mesh/element.h"

#include <vector>

namespace espreso {

enum class CubeFace {
    X_0,
    X_1,
    Y_0,
    Y_1,
    Z_0,
    Z_1,
    NONE
};

enum class CubeEdge {
    X_0_Y_0,
    X_0_Y_1,
    X_1_Y_0,
    X_1_Y_1,
    X_0_Z_0,
    X_0_Z_1,
    X_1_Z_0,
    X_1_Z_1,
    Y_0_Z_0,
    Y_0_Z_1,
    Y_1_Z_0,
    Y_1_Z_1,
    NONE
};

enum class CODE: int;

struct ElementGenerator {

    size_t subelements;
    size_t subnodes[3];
    size_t enodes;
    Element::CODE code;

    virtual void pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const =0;
    virtual void pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeEdge edge) const =0;
    virtual void pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const =0;
    virtual void pushEdge(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeEdge edge) const =0;
    virtual void pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const =0;

    virtual ~ElementGenerator() {}
protected:
    ElementGenerator()
    : subelements(0), subnodes{ 0, 0, 0 }, enodes(0), code(Element::CODE::SIZE) {}
};


}



#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_ELEMENT_H_ */
