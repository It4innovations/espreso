
#ifndef SRC_FETI_PROJECTOR_DUALGRAPH_H_
#define SRC_FETI_PROJECTOR_DUALGRAPH_H_

#include "feti/feti.h"

#include <vector>
#include <map>

namespace espreso {

struct DualGraph {

    void pushVertex(int index, int size);
    void initVertices();

    void setFromDomains(const DecompositionFETI *decomposition, const std::vector<int> &lMap);
    void setFromClusters(const DecompositionFETI *decomposition, const std::vector<int> &lMap);

    struct VertexInfo {
        struct {
            int offset, size;
        } kernel;

        struct {
            struct Indices { int offset, size; };
            std::vector<Indices> indices;
            int total;
        } lambdas;
    //
    //    VertexInfo(): domain(0), koffset(0), kernels(0), ncols(0) {}
    //    VertexInfo(int domain, int koffset, int kernels): domain(domain), koffset(koffset), kernels(kernels), ncols(0) {}
    //
    //    bool operator< (const DomainInfo &other) const { return domain <  other.domain; }
    //    bool operator<=(const DomainInfo &other) const { return domain <= other.domain; }
    //    bool operator!=(const DomainInfo &other) const { return domain != other.domain; }
    };

    std::map<int, VertexInfo> vertices, local;
    std::map<int, std::vector<int> > edges;
};

}

#endif /* SRC_FETI_PROJECTOR_DUALGRAPH_H_ */
