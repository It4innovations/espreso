
#include "logging.h"

namespace espreso {

template<typename Tvalue>
static void gather(const Tvalue &value, Tvalue &min, Tvalue &max, Tvalue &total)
{
	typename std::vector<Tvalue> values(config::MPIsize);

	MPI_Gather(&value, sizeof(Tvalue), MPI_BYTE, values.data(), sizeof(Tvalue), MPI_BYTE, 0, MPI_COMM_WORLD);

	min = max = value;
	total = 0;
	for (size_t i = 0; i < values.size(); i++) {
		total += values[i];
		min = std::min(min, values[i]);
		max = std::max(max, values[i]);
	}
}


template<typename Tvalue>
std::string Info::sumValue(const Tvalue &value)
{
	Tvalue min, max, total;
	gather(value, min, max, total);

	std::stringstream ss;
	ss << total << ", average: " << (double)total / config::MPIsize << " (from " << min << " to " << max << ")";
	return ss.str();
}


template<typename Tvalue>
std::string Info::averageValue(const Tvalue &value)
{
	Tvalue min, max, total;
	gather(value, min, max, total);

	std::stringstream ss;
	ss << "average: " << (double)total / config::MPIsize << " (from " << min << " to " << max << ")";
	return ss.str();
}

}
