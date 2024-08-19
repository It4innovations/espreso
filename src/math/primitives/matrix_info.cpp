
#include "matrix_info.h"

namespace espreso {

#ifdef HAVE_MKL
const int Indexing::CSR = 0;
// const int Indexing::CSR = 1;
#else
const int Indexing::CSR = 0;
#endif

const int Indexing::IJV = 1;

}

