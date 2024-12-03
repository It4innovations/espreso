
#ifndef SRC_FETI_PROJECTOR_DUALGRAPH_H_
#define SRC_FETI_PROJECTOR_DUALGRAPH_H_

#include "feti/feti.h"

#include <vector>
#include <map>

namespace espreso {

struct DualGraph {

    void clear();
    void pushVertex(int index, int size);
    void initVertices();

    void setFromDomains(const DecompositionFETI *decomposition, const std::vector<int> &lMap);
    void setFromClusters(const DecompositionFETI *decomposition, const std::vector<int> &lMap);

    void spread(const DecompositionFETI *decomposition);

    void print();

    struct VertexInfo {
        int rank;

        struct {
            int goffset, loffset, size;
        } kernel;

        struct {
            struct Indices { int offset, size; std::vector<int> lower_neighs; };
            std::vector<Indices> indices;
            int total;
        } lambdas;
    };

    std::map<int, VertexInfo> vertices;
    std::map<int, std::vector<int> > edges;
};

}

#endif /* SRC_FETI_PROJECTOR_DUALGRAPH_H_ */
