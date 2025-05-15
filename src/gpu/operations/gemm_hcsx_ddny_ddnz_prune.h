
#ifndef SRC_GPU_OPERATIONS_GEMM_HCSX_DDNY_DDNZ_PRUNE_H
#define SRC_GPU_OPERATIONS_GEMM_HCSX_DDNY_DDNZ_PRUNE_H

#include "math/primitives_new/allocator_new.h"
#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/vector_dense_data_new.h"
#include "math/operations/prune_csx_matx.h"
#include "gpu/operations/submatrix_ddnx_ddnx_noncontig.h"
#include "gpu/operations/supermatrix_ddnx_ddnx_noncontig.h"
#include "gpu/operations/gemm_dcsx_ddny_ddnz.h"
#include "gpu/operations/gemm_ddnx_ddny_ddnz.h"
#include "gpu/operations/convert_dcsx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class gemm_hcsx_ddny_ddnz_prune
{
public:
    gemm_hcsx_ddny_ddnz_prune() = default;
    gemm_hcsx_ddny_ddnz_prune(const gemm_hcsx_ddny_ddnz_prune &) = delete;
    gemm_hcsx_ddny_ddnz_prune(gemm_hcsx_ddny_ddnz_prune &&) = default;
    gemm_hcsx_ddny_ddnz_prune & operator=(const gemm_hcsx_ddny_ddnz_prune &) = delete;
    gemm_hcsx_ddny_ddnz_prune & operator=(gemm_hcsx_ddny_ddnz_prune &&) = default;
    ~gemm_hcsx_ddny_ddnz_prune() = default;
public:
    void set_config(char spdn_A_, bool prune_rows_, bool prune_cols_);
    void set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_, gpu::dnblas::handle handle_dnblas_);
    void set_matrix_h_A(MatrixCsxView_new<T,I> * h_A_);
    void set_matrix_d_B(MatrixDenseView_new<T> * d_B_);
    void set_matrix_d_C(MatrixDenseView_new<T> * d_C_);
    void set_coefficients(T alpha_, T beta_);
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
    MatrixDenseView_new<T> * d_B = nullptr;
    MatrixDenseView_new<T> * d_C = nullptr;
    T alpha = T{1};
    T beta = T{0};
    char spdn_A = '_';
    bool prune_rows = false;
    bool prune_cols = false;
    void * ws_persistent = nullptr;
    size_t wss_internal = 0;
    size_t wss_persistent = 0;
    size_t wss_tmp_preprocess = 0;
    size_t wss_tmp_perform = 0;
    bool called_set_config = false;
    bool called_set_handles = false;
    bool called_setup = false;
    bool called_preprocess = false;
private:
    std::unique_ptr<AllocatorArena_new> ator_persistent;
    std::unique_ptr<AllocatorArena_new> ator_tmp_linear;
    std::unique_ptr<AllocatorSinglePointer_new> ator_tmp_overlap;
    size_t wss_tmp_preprocess_linear = 0;
    size_t wss_tmp_preprocess_overlap = 0;
    size_t wss_tmp_perform_linear = 0;
    size_t wss_tmp_perform_overlap = 0;
    size_t m = 0;
    size_t n = 0;
    size_t k = 0;
    math::operations::prune_csx_matx<T,I> op_h_prune_A;
    std::unique_ptr<convert_dcsx_ddny<T,I>> op_d_sp2dn_A;
    std::unique_ptr<submatrix_ddnx_ddnx_noncontig<T,I>> op_d_sub_B;
    std::unique_ptr<submatrix_ddnx_ddnx_noncontig<T,I>> op_d_sub_C;
    std::unique_ptr<supermatrix_ddnx_ddnx_noncontig<T,I>> op_d_super_C;
    std::unique_ptr<gemm_dcsx_ddny_ddnz<T,I>> op_gemm_sp;
    std::unique_ptr<gemm_ddnx_ddny_ddnz<T>> op_gemm_dn;
    VectorDenseData_new<I> h_pruned_rows;
    VectorDenseData_new<I> h_pruned_cols;
    VectorDenseData_new<I> d_pruned_rows;
    VectorDenseData_new<I> d_pruned_cols;
    MatrixCsxData_new<T,I> h_A_pruned;
    MatrixCsxData_new<T,I> d_A_pruned_sp;
    MatrixDenseData_new<T> d_A_pruned_dn;
    MatrixDenseData_new<T> d_B_pruned;
    MatrixDenseData_new<T> d_C_pruned;
    MatrixCsxView_new<T,I> * h_A_to_use = nullptr;
    MatrixDenseView_new<T> * d_B_to_use = nullptr;
    MatrixDenseView_new<T> * d_C_to_use = nullptr;

};



}
}
}

#endif /* SRC_GPU_OPERATIONS_GEMM_HCSX_DDNY_DDNZ_PRUNE_H */
