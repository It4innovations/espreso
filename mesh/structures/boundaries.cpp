#include "boundaries.h"

using namespace espreso;

std::ostream& espreso::operator<<(std::ostream& os, const Boundaries &b)
{
	for (size_t i = 0; i < b._boundaries.size(); i++) {
		os << i << ": ";
		for (auto it = b._boundaries[i].begin(); it != b._boundaries[i].end(); ++it) {
			os << *it << " ";
		}
		os << "\n";
	}
	return os;
}




