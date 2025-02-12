
#include "hexahedron20.h"

using namespace espreso;

Hexahedron20Generator::Hexahedron20Generator()
{
    subelements = 1;
    enodes = 20;
    code = Element::CODE::HEXA20;
}

void Hexahedron20Generator::pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const
{
    elements.push_back(indices[ 2]);
    elements.push_back(indices[ 8]);
    elements.push_back(indices[ 6]);
    elements.push_back(indices[ 0]);
    elements.push_back(indices[20]);
    elements.push_back(indices[26]);
    elements.push_back(indices[24]);
    elements.push_back(indices[18]);

    elements.push_back(indices[ 5]);
    elements.push_back(indices[ 7]);
    elements.push_back(indices[ 3]);
    elements.push_back(indices[ 1]);
    elements.push_back(indices[23]);
    elements.push_back(indices[25]);
    elements.push_back(indices[21]);
    elements.push_back(indices[19]);
    elements.push_back(indices[11]);
    elements.push_back(indices[17]);
    elements.push_back(indices[15]);
    elements.push_back(indices[ 9]);
}

void Hexahedron20Generator::pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const
{
    pushNodes(elements, indices, face);
    esize.push_back(8);
    etype.push_back((int)Element::CODE::SQUARE8);
}



