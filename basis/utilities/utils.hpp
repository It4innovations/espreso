

#include "utils.h"


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
	int threads;
	char *var = getenv(name.c_str());
	if (var != NULL) {
		std::stringstream ss(var);
		ss >> value;
	} else {
		ESINFO(eslog::ERROR) << "Set environment variable " << name;
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
		ESINFO(eslog::ERROR) << "Distribution of interval <" << start << "," << end << "> is not possible.";
	}
	size_t size = end - start;
	std::vector<Ttype> distribution(parts + 1, 0);
	size_t chunkSize = std::ceil(size / (double)parts);
	for (size_t t = 1; t < parts; t++) {
		distribution[t] = t * chunkSize - start;
		if (distribution[t] - start > size) {
			distribution[t] = end;
		}
	}
	distribution[parts] = end;

	return distribution;
}
