
#ifndef SRC_FETI_COMMON_DUAL_MAP_H_
#define SRC_FETI_COMMON_DUAL_MAP_H_

#include <vector>

namespace espreso {

template<typename T> struct FETI;

struct Dual_Map {

    template <typename T>
    static void set(FETI<T> &feti);

    struct interval { int start, end, offset, size, total; };

    static int size, total;
    static std::vector<interval> local_intervals;
    static std::vector<int> nmap, neighbors, nsize;
};

}

#endif /* SRC_FETI_COMMON_DUAL_MAP_H_ */
