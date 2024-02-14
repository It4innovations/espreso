
#include "domainsurfacestore.h"
#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"
#include "esinfo/envinfo.h"

#include "mesh/mesh.h"

using namespace espreso;

DomainSurfaceStore::DomainSurfaceStore()
: nodes(NULL),
  enodes(NULL),
  epointers(NULL)
{

}

DomainSurfaceStore::~DomainSurfaceStore()
{
    if (nodes != NULL) { delete nodes; }
    if (enodes != NULL) { delete enodes; }
    if (epointers != NULL) { delete epointers; }
}

size_t DomainSurfaceStore::packedFullSize() const
{
    size_t packedSize = 0;

    packedSize += utils::packedSize(nodes);
    packedSize += utils::packedSize(enodes);
    packedSize += utils::packedSize(dnodes);
    packedSize += utils::packedSize(denodes);
    packedSize += utils::packedSize(coordinates);

    packedSize += utils::packedSize(edistribution);

    packedSize += 1;
    if (epointers != NULL) {
        packedSize += sizeof(size_t) + epointers->datatarray().size() * sizeof(int);
    }

    return packedSize;
}

void DomainSurfaceStore::packFull(char* &p) const
{
    utils::pack(nodes, p);
    utils::pack(enodes, p);
    utils::pack(dnodes, p);
    utils::pack(denodes, p);
    utils::pack(coordinates, p);

    utils::pack(edistribution, p);

    utils::pack(epointers != NULL, p);
    if (epointers != NULL) {
        std::vector<int> eindices;
        eindices.reserve(epointers->datatarray().size());
        for (size_t i = 0; i < epointers->datatarray().size(); ++i) {
            eindices.push_back(static_cast<int>(epointers->datatarray()[i]->code));
        }
        utils::pack(eindices, p);
    }
}

void DomainSurfaceStore::unpackFull(const char* &p)
{
    utils::unpack(nodes, p);
    utils::unpack(enodes, p);
    utils::unpack(dnodes, p);
    utils::unpack(denodes, p);
    utils::unpack(coordinates, p);

    utils::unpack(edistribution, p);

    bool notnull;
    utils::unpack(notnull, p);
    if (notnull) {
        std::vector<int> eindices;
        utils::unpack(eindices, p);
        if (epointers != NULL) {
            delete epointers;
        }
        epointers = new serializededata<esint, Element*>(1, tarray<Element*>(info::env::OMP_NUM_THREADS, eindices.size()));
        for (size_t i = 0; i < eindices.size(); ++i) {
            epointers->datatarray()[i] = &Mesh::edata[eindices[i]];
        }
    }
}
