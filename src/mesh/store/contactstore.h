
#ifndef SRC_MESH_STORE_CONTACTSTORE_H_
#define SRC_MESH_STORE_CONTACTSTORE_H_

#include <cstddef>
#include <vector>

#include "surfacestore.h"

namespace espreso {

#define MIN_SLAVE_COVER_RATIO 0.001

struct Point2D {
    double x, y;

    Point2D(): x(0), y(0) {}
    Point2D(double x, double y): x(x), y(y) {}
    Point2D(const Point &p): x(p.x), y(p.y) {}

    const double* data() const { return &x; }
    double* data() { return &x; }
};

struct Triangle {
    Point p[3];

    Triangle() {}
    Triangle(std::vector<Point> &p, esint p1, esint p2, esint p3): p{ p[p1], p[p2], p[p3] } { }
    Triangle(const Point &p1, const Point &p2, const Point &p3): p{ p1, p2, p3 } { }

    void rotate(const Point &axis, const double &cos, const double &sin)
    {
        p[0].rodrigues(axis, cos, sin);
        p[1].rodrigues(axis, cos, sin);
        p[2].rodrigues(axis, cos, sin);
    }

    double area() const {
        return .5 * ((p[1].x - p[0].x) * (p[2].y - p[0].y) - (p[2].x - p[0].x) * (p[1].y - p[0].y));
    }

    static double area(const Point2D p[]) {
        return .5 * ((p[1].x - p[0].x) * (p[2].y - p[0].y) - (p[2].x - p[0].x) * (p[1].y - p[0].y));
    }
};

struct Interface {
    struct Side {
        esint body, faces, triangleOffset, triangleSize, triangleTotalSize;
        double area;

        Side(esint body): body(body), faces(0), triangleOffset(0), triangleSize(0), triangleTotalSize(0), area(0) {}
    };

    Side from, to;

    Interface(esint from, esint to): from(from), to(to) {}

    void setOrientation()
    {
        if (
            (to.area < from.area && to.faces < from.faces) || // both are smaller
            (1.1 * to.faces < from.faces) || // if the number of faces on the second interface is significantly smaller
            (1.1 * to.area < from.area) // if the number of faces is similar
            ) {

            if (to.faces && to.area) {
                std::swap(from, to);
            }
        }
    }
};

struct SparseSegment {
    esint body;
    esint element;
    esint coordinateOffset;
    esint intersectionOffset;
    esint denseSegmentBegin;
    esint denseSegmentEnd;

    SparseSegment()
    : body(0), element(0), coordinateOffset(0), intersectionOffset(0), denseSegmentBegin(0), denseSegmentEnd(0) {}
    SparseSegment(esint body, esint e, esint coffset, esint ioffset, esint doffset)
    : body(body), element(e), coordinateOffset(coffset), intersectionOffset(ioffset), denseSegmentBegin(doffset), denseSegmentEnd(doffset) {}
};

struct DenseSegment {
    esint neigh;
    esint body;
    esint element;
    esint coordinateOffset;
    esint triangles;
    esint triangleOffset;
    esint skip;

    DenseSegment()
    : neigh(0), body(0), element(0), coordinateOffset(0), triangles(0), triangleOffset(0), skip(false) {}
    DenseSegment(esint n, esint body, esint e, esint coffset, esint toffset)
    : neigh(n), body(body), element(e), coordinateOffset(coffset), triangles(0), triangleOffset(toffset), skip(false) {}
};

struct ContactStore {
    std::vector<int> neighbors, neighborsWithMe;
    std::vector<SurfaceStore*> surfaces; // the last surface is the local surface

    static NodeData *nodeNormals, *nodeMultiplicity;

    serializededata<esint, esint> *pairs;

    serializededata<esint, Triangle>* intersections;

    serializededata<esint, SparseSegment>* sparseSide;
    serializededata<esint, DenseSegment>* denseSide;
    serializededata<esint, Point2D>* planeCoordinates;

    std::vector<Interface> interfaces;

    ContactStore();
    ~ContactStore();

    size_t packedFullSize() const;
    void packFull(char* &p) const;
    void unpackFull(const char* &p);
};

}


#endif /* SRC_MESH_STORE_CONTACTSTORE_H_ */
