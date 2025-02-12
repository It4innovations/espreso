
#include "pyramid13.h"

using namespace espreso;

Pyramid13Generator::Pyramid13Generator()
{
    subelements = 6;
    subnodes[0] = 5;
    subnodes[1] = 5;
    subnodes[2] = 5;
    enodes = 13;
    code = Element::CODE::PYRAMID13;
}

void Pyramid13Generator::pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const
{
    elements.push_back(indices[100]);
    elements.push_back(indices[104]);
    elements.push_back(indices[  4]);
    elements.push_back(indices[  0]);
    elements.push_back(indices[ 62]);

    elements.push_back(indices[102]);
    elements.push_back(indices[ 54]);
    elements.push_back(indices[  2]);
    elements.push_back(indices[ 50]);
    elements.push_back(indices[ 81]);
    elements.push_back(indices[ 83]);
    elements.push_back(indices[ 33]);
    elements.push_back(indices[ 31]);


    elements.push_back(indices[104]);
    elements.push_back(indices[124]);
    elements.push_back(indices[ 24]);
    elements.push_back(indices[  4]);
    elements.push_back(indices[ 62]);

    elements.push_back(indices[114]);
    elements.push_back(indices[ 74]);
    elements.push_back(indices[ 14]);
    elements.push_back(indices[ 54]);
    elements.push_back(indices[ 83]);
    elements.push_back(indices[ 93]);
    elements.push_back(indices[ 43]);
    elements.push_back(indices[ 33]);


    elements.push_back(indices[124]);
    elements.push_back(indices[120]);
    elements.push_back(indices[ 20]);
    elements.push_back(indices[ 24]);
    elements.push_back(indices[ 62]);

    elements.push_back(indices[122]);
    elements.push_back(indices[ 70]);
    elements.push_back(indices[ 22]);
    elements.push_back(indices[ 74]);
    elements.push_back(indices[ 93]);
    elements.push_back(indices[ 91]);
    elements.push_back(indices[ 41]);
    elements.push_back(indices[ 43]);


    elements.push_back(indices[120]);
    elements.push_back(indices[100]);
    elements.push_back(indices[  0]);
    elements.push_back(indices[ 20]);
    elements.push_back(indices[ 62]);

    elements.push_back(indices[110]);
    elements.push_back(indices[ 50]);
    elements.push_back(indices[ 10]);
    elements.push_back(indices[ 70]);
    elements.push_back(indices[ 91]);
    elements.push_back(indices[ 81]);
    elements.push_back(indices[ 31]);
    elements.push_back(indices[ 41]);


    elements.push_back(indices[100]);
    elements.push_back(indices[120]);
    elements.push_back(indices[124]);
    elements.push_back(indices[104]);
    elements.push_back(indices[ 62]);

    elements.push_back(indices[110]);
    elements.push_back(indices[122]);
    elements.push_back(indices[114]);
    elements.push_back(indices[102]);
    elements.push_back(indices[ 81]);
    elements.push_back(indices[ 91]);
    elements.push_back(indices[ 93]);
    elements.push_back(indices[ 83]);


    elements.push_back(indices[  4]);
    elements.push_back(indices[ 24]);
    elements.push_back(indices[ 20]);
    elements.push_back(indices[  0]);
    elements.push_back(indices[ 62]);

    elements.push_back(indices[ 14]);
    elements.push_back(indices[ 22]);
    elements.push_back(indices[ 10]);
    elements.push_back(indices[  2]);
    elements.push_back(indices[ 33]);
    elements.push_back(indices[ 43]);
    elements.push_back(indices[ 41]);
    elements.push_back(indices[ 31]);
}

void Pyramid13Generator::pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeEdge edge) const
{
    return;
}

void Pyramid13Generator::pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const
{
    switch (face) {
    case CubeFace::X_1:
        nodes.push_back(indices[  4]);
        nodes.push_back(indices[ 24]);
        nodes.push_back(indices[124]);
        nodes.push_back(indices[104]);
        nodes.push_back(indices[ 14]);
        nodes.push_back(indices[ 74]);
        nodes.push_back(indices[114]);
        nodes.push_back(indices[ 54]);
        break;
    case CubeFace::Y_1:
        nodes.push_back(indices[ 24]);
        nodes.push_back(indices[ 20]);
        nodes.push_back(indices[120]);
        nodes.push_back(indices[124]);
        nodes.push_back(indices[ 22]);
        nodes.push_back(indices[ 70]);
        nodes.push_back(indices[122]);
        nodes.push_back(indices[ 74]);
        break;
    case CubeFace::X_0:
        nodes.push_back(indices[ 20]);
        nodes.push_back(indices[  0]);
        nodes.push_back(indices[100]);
        nodes.push_back(indices[120]);
        nodes.push_back(indices[ 10]);
        nodes.push_back(indices[ 50]);
        nodes.push_back(indices[110]);
        nodes.push_back(indices[ 70]);
        break;
    case CubeFace::Y_0:
        nodes.push_back(indices[  0]);
        nodes.push_back(indices[  4]);
        nodes.push_back(indices[104]);
        nodes.push_back(indices[100]);
        nodes.push_back(indices[  2]);
        nodes.push_back(indices[ 54]);
        nodes.push_back(indices[102]);
        nodes.push_back(indices[ 50]);
        break;
    case CubeFace::Z_0:
        nodes.push_back(indices[  0]);
        nodes.push_back(indices[ 20]);
        nodes.push_back(indices[ 24]);
        nodes.push_back(indices[  4]);
        nodes.push_back(indices[ 10]);
        nodes.push_back(indices[ 22]);
        nodes.push_back(indices[ 14]);
        nodes.push_back(indices[  2]);
        break;
    case CubeFace::Z_1:
        nodes.push_back(indices[104]);
        nodes.push_back(indices[124]);
        nodes.push_back(indices[120]);
        nodes.push_back(indices[100]);
        nodes.push_back(indices[114]);
        nodes.push_back(indices[122]);
        nodes.push_back(indices[110]);
        nodes.push_back(indices[102]);
        break;
    default:
        break;
    }
}

void Pyramid13Generator::pushEdge(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeEdge edge) const
{
    return;
//    pushNodes(elements, indices, edge);
//    esize.push_back(3);
//    etype.push_back((int)Element::CODE::LINE3);
}

void Pyramid13Generator::pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const
{
    pushNodes(elements, indices, face);
    esize.push_back(8);
    etype.push_back((int)Element::CODE::SQUARE8);
}



