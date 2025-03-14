
#ifndef SRC_GPU_OPERATIONS_TRSM_HCSX_DDNY_DDNY_H
#define SRC_GPU_OPERATIONS_TRSM_HCSX_DDNY_DDNY_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"
#include "gpu/gpu_management.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class trsm_hcsx_ddny_ddny
{
    // solve A * X = B
    // support in-place, where X=B
protected:
    trsm_hcsx_ddny_ddny() = default;
public:
    trsm_hcsx_ddny_ddny(const trsm_hcsx_ddny_ddny &) = delete;
    trsm_hcsx_ddny_ddny(trsm_hcsx_ddny_ddny &&) = delete;
    trsm_hcsx_ddny_ddny & operator=(const trsm_hcsx_ddny_ddny &) = delete;
    trsm_hcsx_ddny_ddny & operator=(trsm_hcsx_ddny_ddny &&) = delete;
    virtual ~trsm_hcsx_ddny_ddny() = default;
public:
    char get_native_place();
    void set_config(char spdn_A);
    void set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_, gpu::dnblas::handle handle_dnblas_);
    void set_matrix_h_A(MatrixCsxView_new<T,I> * h_A_);
    void set_matrix_d_X(MatrixDenseView_new<T> * d_X_);
    void set_matrix_d_B(MatrixDenseView_new<T> * d_B_);
    void setup();
    size_t get_wss_internal();
    size_t get_wss_persistent();
    size_t get_wss_tmp_preprocess();
    size_t get_wss_tmp_perform();
    void set_ws_persistent(void * ws_persistent_);
    void preprocess_submit(void * ws_tmp);
    void perform_submit(void * ws_tmp);
    static void submit_all(gpu::mgm::queue q, gpu::spblas::handle handle_spblas, MatrixCsxView_new<T,I> * h_A, MatrixDenseView_new<T> * d_X, MatrixDenseView_new<T> * d_B, Allocator_new * ator_gpu);
private:
    gpu::mgm::queue q;
    gpu::spblas::handle handle_spblas;
    gpu::dnblas::handle handle_dnblas;
    MatrixCsxView_new<T,I> * h_A = nullptr;
    MatrixDenseView_new<T> * d_X = nullptr;
    MatrixDenseView_new<T> * d_B = nullptr;
    trsm_dcsx_ddny_ddny<T,I> op_d_inner_trsm_sp;
    trsm_ddnx_ddny<T> op_d_inner_trsm_dn;
    convert_dcsx_ddny<T,I> op_d_sp2dn_A;
    MatrixCsxData_new<T,I> d_A_sp;
    MatrixDenseData_new<T,I> d_A_dn;
    void * ws_persistent = nullptr;
    size_t wss_internal = 0; // approximate number of bytes allocated internally
    size_t wss_persistent = 0;
    size_t wss_tmp_preprocess = 0;
    size_t wss_tmp_perform = 0;
    std::unique_ptr<AllocatorArena_new> ator_ws_persistent;
    std::unique_ptr<AllocatorArena_new> ator_ws_tmp_linear;
    std::unique_ptr<AllocatorSinglePointer_new> ator_ws_tmp_overlap;
    size_t wss_tmp_preprocess_linear = 0;
    size_t wss_tmp_preprocess_overlap = 0;
    size_t wss_tmp_perform_linear = 0;
    size_t wss_tmp_perform_overlap = 0;
    char spdn_A = '_';
    bool called_set_config = false;
    bool called_set_handles = false;
    bool called_setup = false;
    bool called_preprocess = false;
};



}
}
}

#endif /* SRC_GPU_OPERATIONS_TRSM_HCSX_DDNY_DDNY_H */
