
#if !defined(HAVE_CUDA) && !defined(HAVE_ROCM)

#include "gpu_dnblas.h"

namespace espreso {
namespace gpu {
namespace dnblas {

    struct _handle {};

    void handle_create(handle & h, mgm::queue & q) {}

    void handle_destroy(handle & h) {}

    void buffer_collect_size(handle & h, size_t & buffersize, const std::function<void(void)> & f) {}

    void buffer_set(handle & h, void * ptr, size_t size) {}

    void buffer_unset(handle & h) {}

    template<typename T, typename I>
    void trsv(handle & h, char mat_symmetry, char transpose, I n, I ld, T * matrix, T * rhs_sol) {}

    template<typename T, typename I>
    void trsm(handle & h, char side, char mat_symmetry, char transpose, I nrows_X, I ncols_X, T * A, I ld_A, T * rhs_sol, I ld_X) {}

    template<typename T, typename I>
    void herk(handle & h, char out_fill, char transpose, I n, I k, T * A, I ld_A, T * C, I ld_C) {}

    template<typename T, typename I>
    void hemv(handle & h, char fill, I n, T * A, I ld_A, T * vec_in, T * vec_out) {}



    #define INSTANTIATE(T,I) \
    template void trsv<T,I>(handle & h, char mat_symmetry, char transpose, I n, I ld, T * matrix, T * rhs_sol); \
    template void trsm<T,I>(handle & h, char side, char mat_symmetry, char transpose, I nrows_X, I ncols_X, T * A, I ld_A, T * rhs_sol, I ld_X); \
    template void herk<T,I>(handle & h, char out_fill, char transpose, I n, I k, T * A, I ld_A, T * C, I ld_C); \
    template void hemv<T,I>(handle & h, char fill, I n, T * A, I ld_A, T * vec_in, T * vec_out);
        // INSTANTIATE(float,                int32_t)
        INSTANTIATE(double,               int32_t)
        // INSTANTIATE(std::complex<float >, int32_t)
        // INSTANTIATE(std::complex<double>, int32_t)
        // INSTANTIATE(float,                int64_t)
        // INSTANTIATE(double,               int64_t)
        // INSTANTIATE(std::complex<float >, int64_t)
        // INSTANTIATE(std::complex<double>, int64_t)
    #undef INSTANTIATE

}
}
}

#endif
