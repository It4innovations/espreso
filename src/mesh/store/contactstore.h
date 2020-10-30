
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

struct ContactStore {

	SurfaceStore *surface;

	std::vector<int> neighbors;
	std::vector<SurfaceStore> halo;

	serializededata<esint, esint> *localPairs;
	serializededata<esint, esint> *neighPairs;

	serializededata<esint, Point>* intersections;

	// [ diagonal, pointer to plane data, nfull, full0[offset,ntria], full1[offset,ntria], ... ]
	serializededata<esint, esint>* interface;
	serializededata<esint, double>* planeData;

	std::vector<ijv> B;

	ContactStore(SurfaceStore *surface);
	~ContactStore();

	size_t packedFullSize() const;
	void packFull(char* &p) const;
	void unpackFull(const char* &p);
};

}


#endif /* SRC_MESH_STORE_CONTACTSTORE_H_ */
