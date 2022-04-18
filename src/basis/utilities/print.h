
#ifndef SRC_BASIS_UTILITIES_PRINT_H_
#define SRC_BASIS_UTILITIES_PRINT_H_

#include "esinfo/mpiinfo.h"
#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "mesh/store/contactstore.h"
#include <iostream>
#include <ostream>
#include <iomanip>
#include <vector>

namespace espreso {

inline std::ostream& operator<<(std::ostream& os, const Point &p) {
//	os << "<" << p.x << " " << p.y << ">";
	os << "<" << p.x << " " << p.y << " " << p.z << ">";
	return os;
}

template<typename T1, typename T2>
std::ostream& operator<< (std::ostream& os, const std::pair<T1, T2> &v)
{
	os << "<" << v.first << ":" << v.second << ">";
	return os;
}

template<typename T, typename A>
std::ostream& operator<< (std::ostream& os, const std::vector<T, A> &v)
{
	for(size_t i = 0; i < v.size(); ++i) {
		os << v[i] << " ";
	}
	os << "\n";
	return os;
}

template <typename TData>
std::ostream& operator<< (std::ostream& os, edata<TData> &data)
{
	os << "[ ";
	for (auto i = data.begin(); i != data.end(); ++i) {
		os << *i << " ";
	}
	os << "]";
	return os;
}

template <typename TEBoundaries, typename TEData>
std::ostream& operator<< (std::ostream& os, const serializededata<TEBoundaries, TEData> &data)
{
	for(auto e = data.cbegin(); e != data.cend(); ++e) {
		os << *e;
	}
	return os;
}

}

#endif /* SRC_BASIS_UTILITIES_PRINT_H_ */
