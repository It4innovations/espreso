
#ifndef SRC_MESH_STORE_DOMAINSURFACESTORE_H_
#define SRC_MESH_STORE_DOMAINSURFACESTORE_H_

#include "basis/containers/point.h"
#include "contactinfo.h"

#include <cstddef>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Element;

struct DomainSurfaceStore {

    // all (for assembler)
    serializededata<esint, esint>* nodes;
    serializededata<esint, esint>* enodes;

    std::vector<size_t> edistribution;

    std::vector<std::vector<int> > dnodes;
    std::vector<std::vector<int> > denodes;
    std::vector<std::vector<Point> > coordinates;

    serializededata<esint, Element*>* epointers;

    DomainSurfaceStore();
    ~DomainSurfaceStore();

    size_t packedFullSize() const;
    void packFull(char* &p) const;
    void unpackFull(const char* &p);
};

}

#endif /* SRC_MESH_STORE_DOMAINSURFACESTORE_H_ */
