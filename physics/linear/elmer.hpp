
#include "linear.h"

namespace physics {

template <>
void Linear<ELMER>::KMf(SparseMatrix &K, SparseMatrix &M, std::vector<double> &f, size_t part, bool dynamics)
{
	std::cerr << "Linear elasticity with ELMER is not implemented yet\n";
	exit(EXIT_FAILURE);
}

}
