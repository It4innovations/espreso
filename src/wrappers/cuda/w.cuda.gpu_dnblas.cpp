
#ifdef HAVE_CUDA

#include "w.cuda.gpu_dnblas.hpp"

namespace espreso {
namespace gpu {
namespace dnblas {

void handle_create(handle & h, mgm::queue & q)
{
    h.inner = new _handle();
    handle_create(*h.inner, *q.inner);
}

void handle_destroy(handle & h)
{
    handle_destroy(*h.inner);
    delete h.inner;
}

void buffer_collect_size(handle & h, size_t & buffersize, const std::function<void(void)> & f)
{
    buffer_collect_size(*h.inner, buffersize, f);
}

void buffer_set(handle & h, void * ptr, size_t size)
{
    buffer_set(*h.inner, ptr, size);
}

void buffer_unset(handle & h)
{
    buffer_unset(*h.inner);
}

template<>
void trsv<double, int>(handle & h, char mat_symmetry, char transpose, int n, int ld, double * matrix, double * rhs_sol)
{
    trsv(*h.inner, mat_symmetry, transpose, n, ld, matrix, rhs_sol);
}

template<>
void trsm<double, int>(handle & h, char side, char mat_symmetry, char transpose, int nrows_X, int ncols_X, double * A, int ld_A, double * rhs_sol, int ld_X)
{
    trsm(*h.inner, side, mat_symmetry, transpose, nrows_X, ncols_X, A, ld_A, rhs_sol, ld_X);
}

template<>
void herk<double, int>(handle & h, char out_fill, char transpose, int n, int k, double * A, int ld_A, double * C, int ld_C)
{
    herk(*h.inner, out_fill, transpose, n, k, A, ld_A, C, ld_C);
}

template<>
void hemv<double, int>(handle & h, char fill, int n, double * A, int ld_A, double * vec_in, double * vec_out)
{
    hemv(*h.inner, fill, n, A, ld_A, vec_in, vec_out);
}

}
}
}

#endif

