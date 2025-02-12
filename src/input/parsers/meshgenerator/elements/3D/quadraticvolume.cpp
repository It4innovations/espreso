
#include "quadraticvolume.h"

using namespace espreso;

QuadraticVolumeGenerator::QuadraticVolumeGenerator()
{
    subnodes[0] = 3;
    subnodes[1] = 3;
    subnodes[2] = 3;
}

void QuadraticVolumeGenerator::pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeEdge edge) const
{
//    switch (edge) {
//    case CubeEdge::X_0_Y_0:
//        nodes.push_back(indices[0]);
//        nodes.push_back(indices[4]);
//        break;
//    case CubeEdge::X_0_Y_1:
//        nodes.push_back(indices[2]);
//        nodes.push_back(indices[6]);
//        break;
//    case CubeEdge::X_0_Z_0:
//        nodes.push_back(indices[2]);
//        nodes.push_back(indices[6]);
//        break;
//    case CubeEdge::X_0_Z_1:
//        nodes.push_back(indices[0]);
//        nodes.push_back(indices[4]);
//        break;
//    case CubeEdge::X_1_Y_0:
//        nodes.push_back(indices[0]);
//        nodes.push_back(indices[1]);
//        break;
//    case CubeEdge::X_1_Y_1:
//        nodes.push_back(indices[5]);
//        nodes.push_back(indices[4]);
//        break;
//    case CubeEdge::X_1_Z_0:
//        nodes.push_back(indices[1]);
//        nodes.push_back(indices[3]);
//        break;
//    case CubeEdge::X_1_Z_1:
//        nodes.push_back(indices[3]);
//        nodes.push_back(indices[7]);
//        break;
//    case CubeEdge::Y_0_Z_0:
//        nodes.push_back(indices[2]);
//        nodes.push_back(indices[6]);
//        break;
//    case CubeEdge::Y_0_Z_1:
//        nodes.push_back(indices[0]);
//        nodes.push_back(indices[4]);
//        break;
//    case CubeEdge::Y_1_Z_0:
//        nodes.push_back(indices[0]);
//        nodes.push_back(indices[1]);
//        break;
//    case CubeEdge::Y_1_Z_1:
//        nodes.push_back(indices[5]);
//        nodes.push_back(indices[4]);
//        break;
//    default:
//        return;
//    }
}

void QuadraticVolumeGenerator::pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const
{

}

void QuadraticVolumeGenerator::pushEdge(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeEdge edge) const
{
    return;
//    pushNodes(elements, indices, edge);
//    esize.push_back(3);
//    etype.push_back((int)Element::CODE::LINE3);
}

void QuadraticVolumeGenerator::pushSquareNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const
{
    switch (face) {
    case CubeFace::X_1:
        nodes.push_back(indices[ 2]);
        nodes.push_back(indices[ 8]);
        nodes.push_back(indices[26]);
        nodes.push_back(indices[20]);

        nodes.push_back(indices[ 5]);
        nodes.push_back(indices[17]);
        nodes.push_back(indices[23]);
        nodes.push_back(indices[11]);
        break;
    case CubeFace::Y_1:
        nodes.push_back(indices[ 8]);
        nodes.push_back(indices[ 6]);
        nodes.push_back(indices[24]);
        nodes.push_back(indices[26]);

        nodes.push_back(indices[ 7]);
        nodes.push_back(indices[15]);
        nodes.push_back(indices[25]);
        nodes.push_back(indices[17]);
        break;
    case CubeFace::X_0:
        nodes.push_back(indices[ 6]);
        nodes.push_back(indices[ 0]);
        nodes.push_back(indices[18]);
        nodes.push_back(indices[24]);

        nodes.push_back(indices[ 3]);
        nodes.push_back(indices[ 9]);
        nodes.push_back(indices[21]);
        nodes.push_back(indices[15]);
        break;
    case CubeFace::Y_0:
        nodes.push_back(indices[ 0]);
        nodes.push_back(indices[ 2]);
        nodes.push_back(indices[20]);
        nodes.push_back(indices[18]);

        nodes.push_back(indices[ 1]);
        nodes.push_back(indices[11]);
        nodes.push_back(indices[19]);
        nodes.push_back(indices[ 9]);
        break;
    case CubeFace::Z_0:
        nodes.push_back(indices[ 0]);
        nodes.push_back(indices[ 6]);
        nodes.push_back(indices[ 8]);
        nodes.push_back(indices[ 2]);

        nodes.push_back(indices[ 3]);
        nodes.push_back(indices[ 7]);
        nodes.push_back(indices[ 5]);
        nodes.push_back(indices[ 1]);
        break;
    case CubeFace::Z_1:
        nodes.push_back(indices[20]);
        nodes.push_back(indices[26]);
        nodes.push_back(indices[24]);
        nodes.push_back(indices[18]);

        nodes.push_back(indices[23]);
        nodes.push_back(indices[25]);
        nodes.push_back(indices[21]);
        nodes.push_back(indices[19]);
        break;
    default:
        return;
    }
}

void QuadraticVolumeGenerator::pushTriangleNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const
{
    switch (face) {
    case CubeFace::X_0:
        nodes.push_back(indices[ 0]);
        nodes.push_back(indices[ 6]);
        nodes.push_back(indices[18]);
        nodes.push_back(indices[ 3]);
        nodes.push_back(indices[12]);
        nodes.push_back(indices[ 9]);

        nodes.push_back(indices[ 6]);
        nodes.push_back(indices[24]);
        nodes.push_back(indices[18]);
        nodes.push_back(indices[15]);
        nodes.push_back(indices[21]);
        nodes.push_back(indices[12]);
        break;
    case CubeFace::X_1:
        nodes.push_back(indices[ 2]);
        nodes.push_back(indices[ 8]);
        nodes.push_back(indices[20]);
        nodes.push_back(indices[ 5]);
        nodes.push_back(indices[14]);
        nodes.push_back(indices[11]);

        nodes.push_back(indices[ 8]);
        nodes.push_back(indices[26]);
        nodes.push_back(indices[20]);
        nodes.push_back(indices[17]);
        nodes.push_back(indices[23]);
        nodes.push_back(indices[14]);
        break;
    case CubeFace::Y_0:
        nodes.push_back(indices[ 0]);
        nodes.push_back(indices[ 2]);
        nodes.push_back(indices[20]);
        nodes.push_back(indices[ 1]);
        nodes.push_back(indices[11]);
        nodes.push_back(indices[10]);

        nodes.push_back(indices[ 0]);
        nodes.push_back(indices[20]);
        nodes.push_back(indices[18]);
        nodes.push_back(indices[10]);
        nodes.push_back(indices[19]);
        nodes.push_back(indices[ 9]);
        break;
    case CubeFace::Y_1:
        nodes.push_back(indices[ 8]);
        nodes.push_back(indices[ 6]);
        nodes.push_back(indices[26]);
        nodes.push_back(indices[ 7]);
        nodes.push_back(indices[16]);
        nodes.push_back(indices[17]);

        nodes.push_back(indices[ 6]);
        nodes.push_back(indices[24]);
        nodes.push_back(indices[26]);
        nodes.push_back(indices[15]);
        nodes.push_back(indices[25]);
        nodes.push_back(indices[16]);
        break;
    case CubeFace::Z_0:
        nodes.push_back(indices[ 0]);
        nodes.push_back(indices[ 6]);
        nodes.push_back(indices[ 2]);
        nodes.push_back(indices[ 3]);
        nodes.push_back(indices[ 4]);
        nodes.push_back(indices[ 1]);

        nodes.push_back(indices[ 6]);
        nodes.push_back(indices[ 8]);
        nodes.push_back(indices[ 2]);
        nodes.push_back(indices[ 7]);
        nodes.push_back(indices[ 5]);
        nodes.push_back(indices[ 4]);
        break;
    case CubeFace::Z_1:
        nodes.push_back(indices[20]);
        nodes.push_back(indices[26]);
        nodes.push_back(indices[24]);
        nodes.push_back(indices[23]);
        nodes.push_back(indices[25]);
        nodes.push_back(indices[22]);

        nodes.push_back(indices[20]);
        nodes.push_back(indices[24]);
        nodes.push_back(indices[18]);
        nodes.push_back(indices[22]);
        nodes.push_back(indices[21]);
        nodes.push_back(indices[19]);
        break;
    default:
        return;
    }
}




