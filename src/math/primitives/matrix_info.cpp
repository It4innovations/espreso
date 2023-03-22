
#include "matrix_info.h"

namespace espreso {

#ifdef USE_SOLVER_SUITESPARSE
const int Indexing::CSR = 0;
#else
const int Indexing::CSR = 1;
#endif

const int Indexing::IJV = 1;

}
