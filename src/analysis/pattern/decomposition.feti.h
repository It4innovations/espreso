
#ifndef SRC_ANALYSIS_PATTERN_DECOMPOSITION_FETI_H_
#define SRC_ANALYSIS_PATTERN_DECOMPOSITION_FETI_H_

#include "decomposition.direct.h"
#include "basis/containers/serializededata.h"

namespace espreso {

//inline bool operator==(const DIndex &left, const DIndex &right) { return left.domain == right.domain && left.index == right.index; }
//inline bool operator!=(const DIndex &left, const DIndex &right) { return !(left == right); }
//inline bool operator <(const DIndex &left, const DIndex &right) { return left.domain == right.domain ? left.index < right.index : left.domain < right.domain; }

struct DecompositionFETI: public DecompositionDirect {
    struct DIndex { int domain, index; };

    int dbegin, dend, dtotal;
    std::vector<int> neighDomain; // first domain index per neighbor, the last is MY OFFSET
    std::vector<int> dsize;

    serializededata<esint, DIndex> *dmap = nullptr;

    bool ismy(int domain) const
    {
        return dbegin <= domain && domain < dend;
    }

    int noffset(int domain) const
    {
        size_t n = 0;
        while (n + 1 < neighbors.size() && neighDomain[n + 1] <= domain) { ++n; }
        return n;
    }

    void update(std::vector<int> &neighbors);

    ~DecompositionFETI()
    {
        if (dmap) delete dmap;
    }
};

}



#endif /* SRC_ANALYSIS_PATTERN_DECOMPOSITION_FETI_H_ */
