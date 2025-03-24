
#ifndef SRC_GPU_OPERATIONS_TRSM_HCSX_DDNY_DDNY_H
#define SRC_GPU_OPERATIONS_TRSM_HCSX_DDNY_DDNY_H

#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_spblas.h"
#include "gpu/gpu_dnblas.h"
#include "gpu/operations/convert_dcsx_ddny.h"
#include "gpu/operations/trsm_dcsx_ddny_ddny.h"
#include "gpu/operations/trsm_ddnx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class trsm_hcsx_ddny_ddny
{
    // solve A * X = B
    // support in-place, where X=B
public:
    trsm_hcsx_ddny_ddny() = default;
    trsm_hcsx_ddny_ddny(const trsm_hcsx_ddny_ddny &) = delete;
    trsm_hcsx_ddny_ddny(trsm_hcsx_ddny_ddny &&) = default;
    trsm_hcsx_ddny_ddny & operator=(const trsm_hcsx_ddny_ddny &) = delete;
    trsm_hcsx_ddny_ddny & operator=(trsm_hcsx_ddny_ddny &&) = default;
    virtual ~trsm_hcsx_ddny_ddny() = default;
public:
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
private:
    gpu::mgm::queue q;
    gpu::spblas::handle handle_spblas;
    gpu::dnblas::handle handle_dnblas;
    MatrixCsxView_new<T,I> * h_A = nullptr;
    MatrixDenseView_new<T> * d_X = nullptr;
    MatrixDenseView_new<T> * d_B = nullptr;
    void * ws_persistent = nullptr;
    size_t wss_internal = 0; // approximate number of bytes allocated internally
    size_t wss_persistent = 0;
    size_t wss_tmp_preprocess = 0;
    size_t wss_tmp_perform = 0;
    char spdn_A = '_';
    bool called_set_config = false;
    bool called_set_handles = false;
    bool called_setup = false;
    bool called_preprocess = false;
private:
    std::unique_ptr<AllocatorArena_new> ator_ws_persistent;
    std::unique_ptr<AllocatorArena_new> ator_ws_tmp_linear;
    std::unique_ptr<AllocatorSinglePointer_new> ator_ws_tmp_overlap;
    size_t wss_tmp_preprocess_linear = 0;
    size_t wss_tmp_preprocess_overlap = 0;
    size_t wss_tmp_perform_linear = 0;
    size_t wss_tmp_perform_overlap = 0;
    MatrixCsxData_new<T,I> d_A_sp;
    MatrixDenseData_new<T> d_A_dn;
    std::unique_ptr<trsm_dcsx_ddny_ddny<T,I>> op_d_inner_trsm_sp;
    std::unique_ptr<trsm_ddnx_ddny<T>> op_d_inner_trsm_dn;
    std::unique_ptr<convert_dcsx_ddny<T,I>> op_d_sp2dn_A;
};



}
}
}

#endif /* SRC_GPU_OPERATIONS_TRSM_HCSX_DDNY_DDNY_H */
