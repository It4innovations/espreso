
#include "esconfig.h"

namespace esconfig {

int MPIrank = 0;
int MPIsize = 1;
Discretization discretization = FEM;

namespace solver {
	double epsilon = .0001;
	size_t maxIterations = 500;
}


}


