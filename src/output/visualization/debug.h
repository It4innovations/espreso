
#ifndef SRC_OUTPUT_VISUALIZATION_SEPARATED_DEBUG_H_
#define SRC_OUTPUT_VISUALIZATION_SEPARATED_DEBUG_H_

#include "writer/vtkwritter.h"
#include "basis/containers/point.h"
#include <string>
#include <sstream>

namespace espreso {

class Mesh;

class DebugOutput {
    static int iteration;

public:
    static void data(
            const std::string &name,
            const std::vector<Point> &points,
            const std::vector<std::vector<esint> > &cells,
            const std::vector<esint> &celltypes,
            const std::vector<std::vector<double> > &celldata);

    static void mesh(double clusterShrinkRatio = .9, double domainShrinkRatio = .95);
    static void faceNeighbors();
    static void meshDual(std::vector<esint> &frames, std::vector<esint> &neighbors);
    static void corners(double clusterShrinkRatio = .9, double domainShrinkRatio = .95);
    static void innerFixPoints(double clusterShrinkRatio = .9, double domainShrinkRatio = .95);
    static void surfaceFixPoints(double clusterShrinkRatio = .9, double domainShrinkRatio = .95);
    static void closeElements(double clusterShrinkRatio = .9, double domainShrinkRatio = .95);
    static void contact(double clusterShrinkRatio = .9, double domainShrinkRatio = .95);
    static void surface(const char* name, double clusterShrinkRatio = .9, double domainShrinkRatio = .95);
    static void warpedNormals(const char* name, const double *displacement = nullptr, double clusterShrinkRatio = .9, double domainShrinkRatio = .95);

protected:
    DebugOutput(double clusterShrinkRatio, double domainShrinkRatio, bool withDomains);
    ~DebugOutput();

    void points(esint nother, esint &noffset, esint &nsize);
    void pointsInDomains(esint nother, esint &noffset, esint &nsize);
    esint elements(esint noffset, esint nother, esint nothernodes);
    esint elementsInDomains(esint noffset, esint nother, esint nothernodes);
    void etypes(esint esize, esint nother);

    std::string _path;
    VTKASCIIWritter _writer;

    Mesh &_mesh;
    Point _ccenter;
    Point *_dcenters;
    double _clusterShrinkRatio, _domainShrinkRatio;
};

}



#endif /* SRC_OUTPUT_VISUALIZATION_SEPARATED_DEBUG_H_ */
