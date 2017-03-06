
#ifndef BASIS_UTILITIES_UTILS_H_
#define BASIS_UTILITIES_UTILS_H_

#include <stdlib.h>
#include <vector>
#include <string>
#include <ostream>

namespace espreso {

struct Esutils
{

	template<typename Ttype>
	static Ttype getEnv(const std::string &name);

	template<typename Ttype>
	static void setFromEnv(Ttype &value, const std::string &name);

	template<typename Ttype>
	static std::vector<Ttype> getDistribution(size_t parts, Ttype size);

	template<typename Ttype>
	static std::vector<Ttype> getDistribution(size_t parts, Ttype start, Ttype end);

	template<typename Ttype>
	static Ttype sizesToOffsets(std::vector<Ttype> &sizes);

	template<typename Ttype>
	static void removeDuplicity(std::vector<Ttype> &elements);

	template<typename Ttype>
	static typename std::vector<Ttype>::const_iterator max_element(const std::vector<Ttype> &elements);

	static std::string createDirectory(const std::vector<std::string> &path);
};


template<typename T>
std::ostream& operator<< (std::ostream& os, const std::vector<T> &v)
{
	for(size_t i = 0; i < v.size(); ++i) {
		os << v[i] << " ";
	}
	os << "\n";
	return os;
}

}

#include "utils.hpp"


#endif /* BASIS_UTILITIES_UTILS_H_ */
