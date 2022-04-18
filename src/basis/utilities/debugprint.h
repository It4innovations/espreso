
#ifndef SRC_BASIS_UTILITIES_DEBUGPRINT_H_
#define SRC_BASIS_UTILITIES_DEBUGPRINT_H_

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include <ostream>
#include <vector>
#include <iomanip>

namespace espreso {

inline std::ostream& operator<<(std::ostream& os, const Point &p) {
	os << "<" << p.x << " " << p.y << " " << p.z << ">\n";
	return os;
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2> &v)
{
	os << "<" << v.first << ":" << v.second << ">\n";
	return os;
}

template<typename T, typename TAlloc>
std::ostream& operator<<(std::ostream& os, const std::vector<T, TAlloc> &v)
{
	os.precision(15);
	os << std::showpos;
	for(size_t i = 0; i < v.size(); ++i) {
		os << std::setw(25) << std::scientific << v[i] << "\n";
	}
	return os;
}

template <typename TData>
std::ostream& operator<<(std::ostream& os, edata<TData> &data)
{
	os << "[ ";
	for (auto i = data.begin(); i != data.end(); ++i) {
		os << *i << " ";
	}
	os << "]\n";
	return os;
}

template <typename TEBoundaries, typename TEData>
std::ostream& operator<<(std::ostream& os, const serializededata<TEBoundaries, TEData> &data)
{
	size_t i = 0;
	for(auto e = data.cbegin(); e != data.cend(); ++e, ++i) {
		os << i << ": " << *e;
	}
	return os;
}

}



#endif /* SRC_BASIS_UTILITIES_DEBUGPRINT_H_ */
