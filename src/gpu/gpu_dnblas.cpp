
#include "gpu_dnblas.h"

namespace espreso {
namespace gpu {
namespace dnblas {

#ifndef HAVE_CUDA

void handle_create(handle & h, mgm::queue & q) {}
void handle_destroy(handle & h) {}

void buffer_collect_size(handle & h, size_t & buffersize, const std::function<void(void)> & f) {}

void buffer_set(handle & h, void * ptr, size_t size) {}

void buffer_unset(handle & h) {}

template<>
void trsv<double, int>(handle & h, char mat_symmetry, char transpose, int n, int ld, double * matrix, double * rhs_sol) {}

template<>
void trsm<double, int>(handle & h, char side, char mat_symmetry, char transpose, int nrows_X, int ncols_X, double * A, int ld_A, double * rhs_sol, int ld_X) {}

template<>
void herk<double, int>(handle & h, char out_fill, char transpose, int n, int k, double * A, int ld_A, double * C, int ld_C) {}

template<>
void hemv<double, int>(handle & h, char fill, int n, double * A, int ld_A, double * vec_in, double * vec_out) {}

#endif

}
}
}
