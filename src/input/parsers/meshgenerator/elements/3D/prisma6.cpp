
#include "prisma6.h"

using namespace espreso;

Prisma6Generator::Prisma6Generator()
{
    subelements = 2;
    enodes = 6;
    code = Element::CODE::PRISMA6;
}

void Prisma6Generator::pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const
{
    elements.push_back(indices[0]);
    elements.push_back(indices[1]);
    elements.push_back(indices[3]);
    elements.push_back(indices[4]);
    elements.push_back(indices[5]);
    elements.push_back(indices[7]);

    elements.push_back(indices[0]);
    elements.push_back(indices[3]);
    elements.push_back(indices[2]);
    elements.push_back(indices[4]);
    elements.push_back(indices[7]);
    elements.push_back(indices[6]);
}

void Prisma6Generator::pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const
{
    switch (face) {
    case CubeFace::X_0:
    case CubeFace::X_1:
    case CubeFace::Y_0:
    case CubeFace::Y_1:
        pushSquareNodes(nodes, indices, face);
        break;
    case CubeFace::Z_0:
    case CubeFace::Z_1:
        pushTriangleNodes(nodes, indices, face);
        break;
    default:
        break;
    }
}

void Prisma6Generator::pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const
{
    pushNodes(elements, indices, face);
    switch (face) {
    case CubeFace::X_0:
    case CubeFace::X_1:
    case CubeFace::Y_0:
    case CubeFace::Y_1:
        esize.push_back(4);
        etype.push_back((int)Element::CODE::SQUARE4);
        break;
    case CubeFace::Z_0:
    case CubeFace::Z_1:
        esize.push_back(3);
        esize.push_back(3);
        etype.push_back((int)Element::CODE::TRIANGLE3);
        etype.push_back((int)Element::CODE::TRIANGLE3);
        break;
    default:
        break;
    }
}
