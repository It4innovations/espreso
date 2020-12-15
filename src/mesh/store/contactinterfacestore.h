
#ifndef SRC_MESH_STORE_CONTACTINTERFACESTORE_H_
#define SRC_MESH_STORE_CONTACTINTERFACESTORE_H_

#include "boundaryregionstore.h"

namespace espreso {

struct ContactInterfaceStore: public BoundaryRegionStore {

	ContactInterfaceStore(const std::string &name): BoundaryRegionStore(name) {}
	ContactInterfaceStore(const char* &packedData): BoundaryRegionStore(packedData) {}
	~ContactInterfaceStore() {};
};

}


#endif /* SRC_MESH_STORE_CONTACTINTERFACESTORE_H_ */
