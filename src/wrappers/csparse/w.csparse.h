#ifndef SRC_WRAPPERS_CSPARSE_W_CSPARSE_H_
#define SRC_WRAPPERS_CSPARSE_W_CSPARSE_H_

#include "feti/generic/SparseMatrix.h"

namespace espreso {
namespace csparse {

// CSparse factorization - 0-based indexing, CSC-format
// Factorize A on CPU, solve mrhs linear system and SpMM on GPU
void CreateLscGpu(SparseMatrix& A, SparseMatrix& B, int order, int tol, int gpu_id, int print_output, SparseMatrix& SC);

}
}

#endif /* SRC_WRAPPERS_CSPARSE_W_CSPARSE_H_ */
