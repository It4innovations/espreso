
#ifndef SRC_FETI_COMMON_DUAL_BUFFER_H_
#define SRC_FETI_COMMON_DUAL_BUFFER_H_

#include <vector>

namespace espreso {

template<typename T> struct FETI;

struct Dual_Map {
    static int nhalo, size;
    static std::vector<int> nmap, neighbors;
};

template <typename T>
struct Dual_Buffer: public Dual_Map {
    template <typename Other>
    static void set(FETI<Other> &feti);
    static std::vector<std::vector<T> > send, recv;
};

template <typename T> std::vector<std::vector<T> > Dual_Buffer<T>::send;
template <typename T> std::vector<std::vector<T> > Dual_Buffer<T>::recv;

}

#endif /* SRC_FETI_COMMON_DUAL_BUFFER_H_ */
