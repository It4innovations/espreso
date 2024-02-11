
#ifndef SRC_FETI_COMMON_DUAL_BUFFER_H_
#define SRC_FETI_COMMON_DUAL_BUFFER_H_

#include <vector>

namespace espreso {

template<typename T> struct FETI;

template <typename T>
struct Dual_Buffer {
    void set(FETI<T> &feti);

	int nhalo, size;
	std::vector<int> nmap, neighbors;
	std::vector<std::vector<T> > send, recv;
};

}

#endif /* SRC_FETI_COMMON_DUAL_BUFFER_H_ */