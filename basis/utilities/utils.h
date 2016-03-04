
#ifndef BASIS_UTILITIES_UTILS_H_
#define BASIS_UTILITIES_UTILS_H_

#include <stdlib.h>
#include "../logging/logging.h"

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
};

#include "utils.hpp"


#endif /* BASIS_UTILITIES_UTILS_H_ */
