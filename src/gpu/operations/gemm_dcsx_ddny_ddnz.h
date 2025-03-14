
#ifndef SRC_GPU_OPERATIONS_GEMM_DCSX_DDNY_DDNZ_H
#define SRC_GPU_OPERATIONS_GEMM_DCSX_DDNY_DDNZ_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_spblas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class gemm_dcsx_ddny_ddnz
{
    // C = alpha * A * B + beta * C
protected:
    gemm_dcsx_ddny_ddnz() = default;
public:
    gemm_dcsx_ddny_ddnz(const gemm_dcsx_ddny_ddnz &) = delete;
    gemm_dcsx_ddny_ddnz(gemm_dcsx_ddny_ddnz &&) = delete;
    gemm_dcsx_ddny_ddnz & operator=(const gemm_dcsx_ddny_ddnz &) = delete;
    gemm_dcsx_ddny_ddnz & operator=(gemm_dcsx_ddny_ddnz &&) = delete;
    virtual ~gemm_dcsx_ddny_ddnz() = 0;
public:
    static std::unique_ptr<gemm_dcsx_ddny_ddnz<T,I>> make();
public:
    void set_handles(gpu::mgm::queue q_, gpu::spblas::handle spblas_handle_);
    void set_matrix_A(MatrixCsxView_new<T,I> * A_);
    void set_matrix_B(MatrixDenseView_new<T> * B_);
    void set_matrix_C(MatrixDenseView_new<T> * C_);
    void set_coefficients(T alpha_, T beta_);
    void setup();
    size_t get_wss_internal();
    size_t get_wss_persistent();
    size_t get_wss_tmp_preprocess();
    size_t get_wss_tmp_perform();
    void set_ws_persistent(void * ws_persistent_);
    void preprocess_submit(void * ws_tmp);
    void perform_submit(void * ws_tmp);
    static void submit_all(gpu::mgm::queue q, gpu::spblas::handle handle_spblas, MatrixCsxView_new<T,I> * A, MatrixDenseView_new<T> * B, MatrixDenseView_new<T> * C, T alpha, T beta, Allocator_new * ator_gpu);
protected:
    gpu::mgm::queue q;
    gpu::spblas::handle handle_spblas;
    MatrixCsxView_new<T,I> * A = nullptr;
    MatrixDenseView_new<T> * B = nullptr;
    MatrixDenseView_new<T> * C = nullptr;
    void * ws_persistent = nullptr;
    size_t wss_internal = 0; // approximate number of bytes allocated internally
    size_t wss_persistent = 0;
    size_t wss_tmp_preprocess = 0;
    size_t wss_tmp_perform = 0;
    T alpha = T{1};
    T beta = T{0};
    bool called_set_handles = false;
    bool called_setup = false;
    bool called_preprocess = false;
protected:
    virtual void internal_setup() {}
    virtual void internal_preprocess(void * /*ws_tmp*/) {}
    virtual void internal_perform(void * /*ws_tmp*/) {}
};

gemm_dcsx_ddny_ddnz::~gemm_dcsx_ddny_ddnz() = default;



}
}
}

#endif /* SRC_GPU_OPERATIONS_GEMM_DCSX_DDNY_DDNZ_H */
