
#ifndef SRC_MESH_STORE_SURFACESTORE_H_
#define SRC_MESH_STORE_SURFACESTORE_H_

#include <cstddef>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Element;

struct SurfaceStore {

	serializededata<esint, esint>* triangles;
	serializededata<esint, esint>* enodes;

	serializededata<esint, esint>* nelements; // only for neighbors definition
	serializededata<esint, esint>* IDs; // only for neighbors definition
	serializededata<esint, esint>* neighbors;

	std::vector<size_t> tdistribution, edistribution;

	esint eoffset;
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