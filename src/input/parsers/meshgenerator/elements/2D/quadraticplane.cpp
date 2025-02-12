
#include "quadraticplane.h"

using namespace espreso;

QuadraticPlaneGenerator::QuadraticPlaneGenerator()
{
    subnodes[0] = 3;
    subnodes[1] = 3;
    subnodes[2] = 1;
}

void QuadraticPlaneGenerator::pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeEdge edge) const
{
    switch (edge) {
    case CubeEdge::X_0_Z_0:
        nodes.push_back(indices[6]);
        nodes.push_back(indices[0]);
        nodes.push_back(indices[3]);
        break;
    case CubeEdge::X_1_Z_0:
        nodes.push_back(indices[2]);
        nodes.push_back(indices[8]);
        nodes.push_back(indices[5]);
        break;
    case CubeEdge::Y_0_Z_0:
        nodes.push_back(indices[0]);
        nodes.push_back(indices[2]);
        nodes.push_back(indices[1]);
        break;
    case CubeEdge::Y_1_Z_0:
        nodes.push_back(indices[8]);
        nodes.push_back(indices[6]);
        nodes.push_back(indices[7]);
        break;
    default:
        return;
    }
}

void QuadraticPlaneGenerator::pushEdge(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeEdge edge) const
{
    pushNodes(elements, indices, edge);
    esize.push_back(3);
    etype.push_back((int)Element::CODE::LINE3);
}

void QuadraticPlaneGenerator::pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const
{
    nodes.push_back(indices[0]);
    nodes.push_back(indices[1]);
    nodes.push_back(indices[2]);
    nodes.push_back(indices[3]);
    nodes.push_back(indices[4]);
    nodes.push_back(indices[5]);
    nodes.push_back(indices[6]);
    nodes.push_back(indices[7]);
    nodes.push_back(indices[8]);
    return;
}

void QuadraticPlaneGenerator::pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const
{
    return;
}

