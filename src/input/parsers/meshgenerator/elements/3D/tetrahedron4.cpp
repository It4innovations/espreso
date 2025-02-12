
#include "tetrahedron4.h"

using namespace espreso;

Tetrahedron4Generator::Tetrahedron4Generator()
{
    subelements = 6;
    enodes = 4;
    code = Element::CODE::TETRA4;
}

void Tetrahedron4Generator::pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const
{
    elements.push_back(indices[0]);
    elements.push_back(indices[3]);
    elements.push_back(indices[2]);
    elements.push_back(indices[4]);

    elements.push_back(indices[3]);
    elements.push_back(indices[2]);
    elements.push_back(indices[4]);
    elements.push_back(indices[6]);

    elements.push_back(indices[7]);
    elements.push_back(indices[3]);
    elements.push_back(indices[4]);
    elements.push_back(indices[6]);

    elements.push_back(indices[3]);
    elements.push_back(indices[5]);
    elements.push_back(indices[7]);
    elements.push_back(indices[4]);

    elements.push_back(indices[1]);
    elements.push_back(indices[5]);
    elements.push_back(indices[3]);
    elements.push_back(indices[4]);

    elements.push_back(indices[0]);
    elements.push_back(indices[4]);
    elements.push_back(indices[1]);
    elements.push_back(indices[3]);
}

void Tetrahedron4Generator::pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const
{
    pushTriangleNodes(elements, indices, face);
    esize.push_back(3);
    esize.push_back(3);
    etype.push_back((int)Element::CODE::TRIANGLE3);
    etype.push_back((int)Element::CODE::TRIANGLE3);
}



