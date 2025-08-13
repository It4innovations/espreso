
#ifndef SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_
#define SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_

#include "regionstore.h"
#include "nameddata.h"
#include "nodestore.h"

#include <cstddef>
#include <vector>
#include <string>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Element;

struct BoundaryElementData: public NamedData {
    enum class Type {
        NODES,
        ELEMENTS
    };

    BoundaryElementData(int dimension, Type type, DataType datatype, const std::string &name): NamedData(dimension, datatype, name) {}
};

struct BoundaryRegionStore: public RegionStore {
    int originalDimension, dimension;
    double area;

    NodeData *nodeNormals, *nodeMultiplicity;

    BoundaryElementData* appendData(int dimension, BoundaryElementData::Type type, NamedData::DataType datatype, const std::string &name = "", step::TYPE restriction = step::TYPE::TIME, bool toOutput = true);

    serializededata<esint, esint>* elements;
    serializededata<esint, esint>* triangles;

    serializededata<esint, Element*>* epointers;
    serializededata<esint, esint>* emembership;

    std::vector<ElementsInterval> eintervals;
    std::vector<esint> eintervalsDistribution;

    std::vector<BoundaryElementData*> data;

    size_t packedFullSize() const;
    void packFull(char* &p) const;
    void unpackFull(const char* &p);

    size_t packedSize() const;
    void pack(char* &p) const;
    void unpack(const char* &p);

    void permute(const std::vector<esint> &permutation) { permute(permutation, distribution.threads); }
    void permute(const std::vector<esint> &permutation, const std::vector<size_t> &threading);

    BoundaryRegionStore(const std::string &name);
    BoundaryRegionStore(const char* &packedData);
    ~BoundaryRegionStore();
};

}


#endif /* SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_ */
