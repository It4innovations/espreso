
#include "../logging/logging.h"
#include <cmath>

#include "utils.h"

namespace espreso {

template<typename Ttype>
Ttype Esutils::getEnv(const std::string &name)
{
	Ttype value;
	setFromEnv(value, name);
	return value;
}

template<typename Ttype>
void Esutils::setFromEnv(Ttype &value, const std::string &name)
{
	char *var = getenv(name.c_str());
	if (var != NULL) {
		std::stringstream ss(var);
		ss >> value;
	} else {
		ESINFO(ERROR) << "Set environment variable " << name;
	}
}


template<typename Ttype>
std::vector<Ttype> Esutils::getDistribution(size_t parts, Ttype size)
{
	return getDistribution(parts, (Ttype)0, size);
}

template<typename Ttype>
std::vector<Ttype> Esutils::getDistribution(size_t parts, Ttype start, Ttype end)
{
	if (start > end) {
		ESINFO(ERROR) << "Distribution of interval <" << start << "," << end << "> is not possible.";
	}
	size_t size = end - start;
	std::vector<Ttype> distribution(parts + 1, 0);
	size_t chunkSize = std::ceil(size / (double)parts);
	for (size_t t = 1; t < parts; t++) {
		distribution[t] = t * chunkSize + start;
		if (distribution[t] > size) {
			distribution[t] = end;
		}
	}
	distribution[0] = start;
	distribution[parts] = end;

	return distribution;
}

template<typename Ttype>
Ttype Esutils::sizesToOffsets(std::vector<Ttype> &sizes)
{
	Ttype sum = 0;
	for (size_t i = 0; i < sizes.size(); i++) {
		Ttype tmp = sizes[i];
		sizes[i] = sum;
		sum += tmp;
	}
	return sum;
}

template<typename Ttype>
void Esutils::removeDuplicity(std::vector<Ttype> &elements)
{
	if (elements.size() == 0) {
		return;
	}
	size_t unique = 0;
	for (size_t d = 1; d < elements.size(); d++) {
		if (elements[unique] != elements[d]) {
			elements[++unique] = elements[d];
		}
	}

	elements.resize(unique + 1);
}

template<typename Ttype>
typename std::vector<Ttype>::const_iterator Esutils::max_element(const std::vector<Ttype> &elements)
{
	auto max = elements.begin();
	for (size_t i = 1; i < elements.size(); i++) {
		if (*max < elements[i]) {
			max = elements.begin() + i;
		}
	}
	return max;
}

}
