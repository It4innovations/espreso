
#ifndef SRC_GPU_OPERATIONS_TRSM_DDNX_DDNY_H
#define SRC_GPU_OPERATIONS_TRSM_DDNX_DDNY_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_dnblas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
class trsm_ddnx_ddny
{
    // solve A * X = B
protected:
    trsm_ddnx_ddny() = default;
public:
    trsm_ddnx_ddny(const trsm_ddnx_ddny &) = delete;
    trsm_ddnx_ddny(trsm_ddnx_ddny &&) = delete;
    trsm_ddnx_ddny & operator=(const trsm_ddnx_ddny &) = delete;
    trsm_ddnx_ddny & operator=(trsm_ddnx_ddny &&) = delete;
    virtual ~trsm_ddnx_ddny() = 0;
public:
    static std::unique_ptr<trsm_ddnx_ddny<T>> make();
public:
    void set_handles(gpu::mgm::queue q_, gpu::dnblas::handle dnblas_handle_);
    void set_matrix_A(MatrixDenseView_new<T,I> A_);
    void set_matrix_X(MatrixDenseView_new<T> X_);
    void setup();
    size_t get_wss_tmp_perform();
    void perform_submit(void * ws_tmp);
    static void submit_all(gpu::mgm::queue q, gpu::dnblas::handle handle_dnblas, MatrixDenseView_new<T> A, MatrixDenseView_new<T> X, Allocator_new * ator_gpu);
protected:
    gpu::mgm::queue q;
    gpu::dnblas::handle handle_dnblas;
    MatrixDenseView_new<T> A;
    MatrixDenseView_new<T> X;
    size_t wss_tmp_perform = 0;
    bool called_set_handles = false;
    bool called_set_A = false;
    bool called_set_X = false;
    bool called_setup = false;
protected:
    virtual void internal_set_matrix_A() {}
    virtual void internal_set_matrix_X() {}
    virtual void internal_setup() {}
    virtual void internal_perform(void * /*ws_tmp*/) {}
};

trsm_ddnx_ddny::~trsm_ddnx_ddny() = default;



}
}
}

#endif /* SRC_GPU_OPERATIONS_TRSM_DDNX_DDNY_H */
