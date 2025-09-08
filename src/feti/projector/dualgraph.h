
#ifndef SRC_FETI_PROJECTOR_DUALGRAPH_H_
#define SRC_FETI_PROJECTOR_DUALGRAPH_H_

#include "feti/feti.h"

#include <vector>
#include <map>

namespace espreso {

struct DualGraph {

    void clear();
    void pushVertex(int offset, int index, int size);

    void set(const DecompositionFETI *decomposition, const std::vector<int> &lMap);
    void spread(const DecompositionFETI *decomposition);

    void print();

    struct DomainVertexInfo {
        int offset; // offset to local K
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

    struct ClusterVertexInfo: DomainVertexInfo {
        std::vector<int> offset;
    };

    struct {
        std::map<int, DomainVertexInfo> vertices;
        std::map<int, std::vector<int> > edges;
    } domains;

    struct {
        std::map<int, ClusterVertexInfo> vertices;
        std::map<int, std::vector<int> > edges;
    } clusters;
};

}

#endif /* SRC_FETI_PROJECTOR_DUALGRAPH_H_ */
