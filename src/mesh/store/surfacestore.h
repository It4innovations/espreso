
#ifndef SRC_MESH_STORE_SURFACESTORE_H_
#define SRC_MESH_STORE_SURFACESTORE_H_

#include "basis/containers/point.h"
#include "contactinfo.h"

#include <cstddef>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Element;

struct SurfaceStore {

	serializededata<esint, esint>* parents;
	serializededata<esint, ContactInfo>* contact;
	serializededata<esint, esint>* body;
	serializededata<esint, esint>* triangles;
	serializededata<esint, esint>* nodes;
	serializededata<esint, esint>* nIDs;
	serializededata<esint, Point>* coordinates;
	serializededata<esint, esint>* enodes;

	// from coarse nodes + warping in the case of squares
	// base + normal -> normal implicit definition
	// base + parameters -> parametric definition of a wrapping polygon
	serializededata<esint, Point>* base;
	serializededata<esint, Point>* normal;
	serializededata<esint, Point>* parameters;

	serializededata<esint, esint>* nelements; // only for neighbors definition
	serializededata<esint, esint>* IDs; // only for neighbors definition
	serializededata<esint, esint>* neighbors;

	std::vector<size_t> tdistribution, edistribution;

	esint offset;
	esint size;
	esint totalSize;
	serializededata<esint, Element*>* epointers;
	std::vector<esint> ecounters;

	SurfaceStore();
	~SurfaceStore();

	size_t packedFullSize() const;
	void packFull(char* &p) const;
	void unpackFull(const char* &p);
};

}


#endif /* SRC_MESH_STORE_SURFACESTORE_H_ */
