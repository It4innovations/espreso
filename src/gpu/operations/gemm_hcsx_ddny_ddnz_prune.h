
#ifndef SRC_GPU_OPERATIONS_GEMM_HCSX_DDNY_DDNZ_PRUNE_H
#define SRC_GPU_OPERATIONS_GEMM_HCSX_DDNY_DDNZ_PRUNE_H

#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/vector_dense_data_new.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class gemm_hcsx_ddny_ddnz_prune
{
public:
    gemm_csx_dny_dny_prune() = default;
    gemm_csx_dny_dny_prune(const gemm_csx_dny_dny_prune &) = delete;
    gemm_csx_dny_dny_prune(gemm_csx_dny_dny_prune &&) = delete;
    gemm_csx_dny_dny_prune & operator=(const gemm_csx_dny_dny_prune &) = delete;
    gemm_csx_dny_dny_prune & operator=(gemm_csx_dny_dny_prune &&) = delete;
    ~gemm_csx_dny_dny_prune();
public:
    void set_config(char spdn_A_, bool prune_rows_, bool prune_cols_);
    void set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_, gpu::dnblas::handle handle_dnblas_);
    void set_matrix_h_A(MatrixCsxView_new<T,I> h_A_);
    void set_matrix_d_B(MatrixDenseView_new<T> d_B_);
    void set_matrix_d_C(MatrixDenseView_new<T> d_C_);
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
    MatrixCsxView_new<T,I> h_A;
    MatrixDenseView_new<T> d_B;
    MatrixDenseView_new<T> d_C;
    T alpha = T{1};
    T beta = T{0};
    char spdn_A = '_';
    bool prune_rows = false;
    bool prune_cols = false;
private:
    void * ws_persistent = nullptr;
    std::unique_ptr<AllocatorArena_new> ator_persistent;
    size_t wss_internal = 0;
    size_t wss_persistent = 0;
    size_t wss_tmp_preprocess = 0;
    size_t wss_tmp_perform = 0;
    size_t m = 0;
    size_t n = 0;
    size_t k = 0;
    prune_csx_matx<T,I> op_h_prune_A;
    VectorDenseData_new<I> pruned_rows;
    VectorDenseData_new<I> pruned_cols;
    MatrixCsxData_new<T,I> h_A_pruned;
    // MatrixCsxData_new<T,I> d_A_pruned_sp;
    // MatrixDenseData_new<T> d_A_pruned_dn;
    // MatrixDenseData_new<T> d_B_pruned;
    // MatrixDenseData_new<T> d_C_pruned;
    MatrixCsxView_new<T,I> * h_A_to_use = nullptr;
    // MatrixDenseView_new<T> * d_B_to_use = nullptr;
    // MatrixDenseView_new<T> * d_C_to_use = nullptr;
    bool called_set_config = false;
    bool called_set_handles = false;
    bool called_set_A = false;
    bool called_set_B = false;
    bool called_set_C = false;
    bool called_setup = false;
    bool called_preprocess = false;

};



}
}
}

#endif /* SRC_GPU_OPERATIONS_GEMM_HCSX_DDNY_DDNZ_PRUNE_H */
