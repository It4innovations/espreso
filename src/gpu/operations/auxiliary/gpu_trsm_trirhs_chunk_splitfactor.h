
#ifndef SRC_GPU_OPERATIONS_AUXILIARY_GPU_TRSM_TRIRHS_CHUNK_SPLITFACTOR_H
#define SRC_GPU_OPERATIONS_AUXILIARY_GPU_TRSM_TRIRHS_CHUNK_SPLITFACTOR_H

#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/vector_dense_view_new.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class gpu_trsm_trirhs_chunk_splitfactor
{
public:
    struct config
    {
        char trsm_factor_spdn = '_'; // Sparse, Dense
        char trsm_factor_order = '_'; // Rowmajor, Colmajor
        char gemm_factor_spdn = '_'; // Sparse, Dense
        char gemm_factor_order = '_'; // Rowmajor, Colmajor
        char gemm_factor_prune = '_'; // No, Rows only, Cols only, All
    };
public:
    gpu_trsm_trirhs_chunk_splitfactor() = default;
    gpu_trsm_trirhs_chunk_splitfactor(const gpu_trsm_trirhs_chunk_splitfactor &) = delete;
    gpu_trsm_trirhs_chunk_splitfactor(gpu_trsm_trirhs_chunk_splitfactor &&) = delete;
    gpu_trsm_trirhs_chunk_splitfactor & operator=(const gpu_trsm_trirhs_chunk_splitfactor &) = delete;
    gpu_trsm_trirhs_chunk_splitfactor & operator=(gpu_trsm_trirhs_chunk_splitfactor &&) = delete;
    ~gpu_trsm_trirhs_chunk_splitfactor() = default;
public:
    void set_config(config cfg_);
    void set_range(size_t k_start_, size_t k_end_);
    void set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_, gpu::dnblas::handle handle_dnblas_);
    void set_matrix_h_L(MatrixCsxView_new<T,I> * h_L_);
    void set_matrix_d_X(MatrixDenseView_new<T> * d_X_);
    void set_h_X_rowtrails(VectorDenseView_new<I> * h_X_rowtrails_);
    void setup();
    size_t get_wss_internal();
    size_t get_wss_persistent();
    size_t get_wss_tmp_preprocess();
    size_t get_wss_tmp_perform();
    void set_ws_persistent(void * ws_persistent_);
    void preprocess_submit(void * ws_tmp);
    void perform_submit(void * ws_tmp);
private:
    config cfg;
    size_t k_start = 0;
    size_t k_end = 0;
    size_t k_size = 0;
    size_t rhs_end = 0;
    gpu::mgm::queue q;
    gpu::spblas::handle handle_spblas;
    gpu::dnblas::handle handle_dnblas;
    MatrixCsxView_new<T,I> * h_L = nullptr;
    MatrixDenseView_new<T> * d_X = nullptr;
    VectorDenseView_new<I> * h_X_rowtrails = nullptr;
    void * ws_persistent = nullptr;
    size_t wss_internal = 0;
    size_t wss_persistent = 0;
    size_t wss_tmp_preprocess = 0;
    size_t wss_tmp_perform = 0;
    bool called_set_config = false;
    bool called_set_range = false;
    bool called_set_handles = false;
    bool called_setup = false;
    bool called_preprocess = false;
private:
    std::unique_ptr<AllocatorArena_new> ator_ws_persistent;
    std::unique_ptr<AllocatorArena_new> ator_ws_tmp_linear;
    std::unique_ptr<AllocatorSinglePointer_new> ator_ws_tmp_overlap;
    size_t wss_tmp_preprocess_linear = 0;
    size_t wss_tmp_preprocess_overlap = 0;
    size_t wss_tmp_peform_linear = 0;
    size_t wss_tmp_peform_overlap = 0;
    math::operations::submatrix_csx_csy_map<T,I> op_h_submatrix_L_top;
    math::operations::submatrix_csx_csy_map<T,I> op_h_submatrix_L_bot;
    math::operations::submatrix_dnx_dnx_view<T> op_submatrix_d_X_top;
    math::operations::submatrix_dnx_dnx_view<T> op_submatrix_d_X_bot;
    trsm_hcsx_ddny_ddny<T,I> op_trsm;
    gemm_hcsx_ddny_ddnz_prune<T,I> op_gemm;
    MatrixCsxData_new<T,I> h_sub_L_top_sp;
    MatrixCsxData_new<T,I> h_sub_L_bot_sp;
    MatrixDenseView_new<T> d_sub_X_top;
    MatrixDenseView_new<T> d_sub_X_bot;
};



}
}
}

#endif /* SRC_GPU_OPERATIONS_AUXILIARY_GPU_TRSM_TRIRHS_CHUNK_SPLITFACTOR_H */
