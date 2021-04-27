
#ifndef SRC_MESH_STORE_REGIONSTORE_H_
#define SRC_MESH_STORE_REGIONSTORE_H_

#include "elementinfo.h"
#include "elementsinterval.h"
#include "nodeuniquenessinfo.h"

#include <vector>
#include <string>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct RegionStore {

	std::string name;

	ElementsDistributionInfo distribution;

	serializededata<esint, esint>* nodes;
	NodeUniquenessInfo nodeInfo;

	esint getPosition(esint node) const;

	size_t packedFullSize() const;
	void packFull(char* &p) const;
	void unpackFull(const char* &p);

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	RegionStore(const std::string &name);
	RegionStore(const char* &packedData);
	~RegionStore();
};

}




#endif /* SRC_MESH_STORE_REGIONSTORE_H_ */

