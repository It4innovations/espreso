
#ifndef SRC_MESH_STORE_ELEMENTSREGIONSTORE_H_
#define SRC_MESH_STORE_ELEMENTSREGIONSTORE_H_

#include "regionstore.h"
#include "contactinfo.h"

#include <vector>
#include <string>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct SurfaceStore;

struct ElementsRegionStore: public RegionStore {
	serializededata<esint, esint>* elements;

	std::vector<ElementsInterval> eintervals;

	SurfaceStore *surface;

	std::vector<int> bodies;
	std::vector<esint> bodyElements;
	std::vector<esint> bodyFaces;
	ContactInfo contact;

	Point orientation; // temporary until new assembler updates

	size_t packedFullSize() const;
	void packFull(char* &p) const;
	void unpackFull(const char* &p);

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	ElementsRegionStore(const std::string &name);
	ElementsRegionStore(const char* &packedData);
	~ElementsRegionStore();
};

}


#endif /* SRC_MESH_STORE_ELEMENTSREGIONSTORE_H_ */
