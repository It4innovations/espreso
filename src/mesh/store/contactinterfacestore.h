
#ifndef SRC_MESH_STORE_CONTACTINTERFACESTORE_H_
#define SRC_MESH_STORE_CONTACTINTERFACESTORE_H_

#include "boundaryregionstore.h"

namespace espreso {

struct ContactInterfaceStore: public BoundaryRegionStore {

    esint interfaceIndex;

    ContactInterfaceStore(const std::string &name, esint interfaceIndex);
    ContactInterfaceStore(const char* &packedData);
    ~ContactInterfaceStore() {};

    size_t packedFullSize() const;
    void packFull(char* &p) const;
    void unpackFull(const char* &p);
};

}


#endif /* SRC_MESH_STORE_CONTACTINTERFACESTORE_H_ */
