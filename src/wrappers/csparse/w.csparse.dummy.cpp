#include "w.csparse.h"
#include "esinfo/eslog.h"

#ifndef HAVE_CSPARSE

namespace espreso {

void csparse::CreateLscGpu(SparseMatrix& A, SparseMatrix& B, int order, int tol, int gpu_id, int print_output, SparseMatrix& SC) {
    eslog::warning("ESPRESO run-time error: cannot call CSparse library (the library is not linked).\n");
}

}

#endif