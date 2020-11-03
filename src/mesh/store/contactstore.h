
#ifndef SRC_MESH_STORE_CONTACTSTORE_H_
#define SRC_MESH_STORE_CONTACTSTORE_H_

#include <cstddef>
#include <vector>

#include "surfacestore.h"

namespace espreso {

struct ijv {
	int i , j;
	double v;
	ijv(): i(0), j(0), v(0) {}
	ijv(int i, int j, double v): i(i), j(j), v(v) {}

	bool operator<(ijv &other) { return i == other.i ? j < other.j : i < other.i; }
	bool operator==(ijv &other) { return i == other.i && j == other.j; }
	bool operator!=(ijv &other) { return !(*this == other); }
};

struct Triangle {
	Point p[3];

	Triangle() {}
	Triangle(std::vector<Point> &p, esint p1, esint p2, esint p3): p{ p[p1], p[p2], p[p3] } { }
	Triangle(const Point &p1, const Point &p2, const Point &p3): p{ p1, p2, p3 } { }

	void rotate(const Point &center, const Point &axis, const double &cos, const double &sin)
	{
		p[0] += center;
		p[1] += center;
		p[2] += center;
		p[0].rodrigues(axis, cos, sin);
		p[1].rodrigues(axis, cos, sin);
		p[2].rodrigues(axis, cos, sin);
	}

	double area() {
		return .5 * ((p[1].x - p[0].x) * (p[2].y - p[0].y) - (p[2].x - p[0].x) * (p[1].y - p[0].y));
	}
};

struct Interface {
	struct Side {
		esint body, faces;
		double area;

		Side(esint body): body(body), faces(0), area(0) {}
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

			std::swap(from, to);
		}
	}
};

struct ContactStore {
	std::vector<int> neighbors, neighborsWithMe;
	std::vector<SurfaceStore*> surfaces; // the last surface is the local surface

	serializededata<esint, esint> *pairs;

	serializededata<esint, Triangle>* intersections;

	// [ diagonal, pointer to plane data, nfull, full0[offset,ntria], full1[offset,ntria], ... ]
	serializededata<esint, esint>* interface;
	serializededata<esint, double>* planeData;

	std::vector<Interface> interfaces;

	std::vector<ijv> B;

	ContactStore();
	~ContactStore();

	size_t packedFullSize() const;
	void packFull(char* &p) const;
	void unpackFull(const char* &p);
};

}


#endif /* SRC_MESH_STORE_CONTACTSTORE_H_ */
