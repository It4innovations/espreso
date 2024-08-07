
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
    void trsv(handle & h, I n, T * A, I ld_A, char order_A, char op_A, char fill_A, T * x) {}

    template<typename T, typename I>
    void trsm(handle & h, char side, I n, I nrhs, T * A, I ld_A, char order_A, char op_A, char fill_A, T * X, I ld_X, char order_X, char op_X) {}

    template<typename T, typename I>
    void herk(handle & h, I n, I k, T * A, I ld_A, char order_A, char op_A, T * C, I ld_C, char order_C, char fill_C) {}

    template<typename T, typename I>
    void hemv(handle & h, I n, T * A, I ld_A, char order_A, char op_A, char fill_A, T * x, T * y) {}

    template<typename T, typename I>
    void gemv(handle & h, I m, I n, T * A, I ld_A, char order_A, char op_A, T * x, T * y) {}

}
}
}

#include "gpu/gpu_dnblas.inst.hpp"

#endif
