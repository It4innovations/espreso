
#ifndef SRC_GPU_OPERATIONS_AUXILIARY_GPU_TRSM_TRIRHS_CHUNK_SPITRHS_H
#define SRC_GPU_OPERATIONS_AUXILIARY_GPU_TRSM_TRIRHS_CHUNK_SPITRHS_H

#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/operations/submatrix_dnx_dnx_view.h"
#include "gpu/operations/submatrix_dcsx_dcsx.h"
#include "gpu/operations/trsm_dcsx_ddny_ddny.h"
#include "gpu/operations/trsm_ddnx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class gpu_trsm_trirhs_chunk_splitrhs
{
public:
    gpu_trsm_trirhs_chunk_splitrhs() = default;
    gpu_trsm_trirhs_chunk_splitrhs(const gpu_trsm_trirhs_chunk_splitrhs &) = delete;
    gpu_trsm_trirhs_chunk_splitrhs(gpu_trsm_trirhs_chunk_splitrhs &&) = default;
    gpu_trsm_trirhs_chunk_splitrhs & operator=(const gpu_trsm_trirhs_chunk_splitrhs &) = delete;
    gpu_trsm_trirhs_chunk_splitrhs & operator=(gpu_trsm_trirhs_chunk_splitrhs &&) = default;
    ~gpu_trsm_trirhs_chunk_splitrhs() = default;
public:
    void set_config(char factor_spdn_);
    void set_range(size_t rhs_start_, size_t rhs_end_);
    void set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_, gpu::dnblas::handle handle_dnblas_);
    void set_matrix_d_L_sp(MatrixCsxView_new<T,I> * d_L_sp_);
    void set_matrix_d_L_dn(MatrixDenseView_new<T> * d_L_dn_);
    void set_matrix_d_X(MatrixDenseView_new<T> * d_X_);
    void set_h_X_colpivots(VectorDenseView_new<I> * h_X_colpivots_);
    void set_h_L_nnzinsubs(VectorDenseView_new<I> * h_L_nnzinsubs_);
    void setup();
    size_t get_wss_internal();
    size_t get_wss_persistent();
    size_t get_wss_tmp_preprocess();
    size_t get_wss_tmp_perform();
    void set_ws_persistent(void * ws_persistent_);
    void preprocess_submit(void * ws_tmp);
    void perform_submit(void * ws_tmp);
private:
    char factor_spdn = '_';
    size_t rhs_start = 0;
    size_t rhs_end = 0;
    size_t rhs_size = 0;
    size_t k_start = 0;
    size_t k_size = 0;
    gpu::mgm::queue q;
    gpu::spblas::handle handle_spblas;
    gpu::dnblas::handle handle_dnblas;
    MatrixCsxView_new<T,I> * d_L_sp = nullptr;
    MatrixDenseView_new<T> * d_L_dn = nullptr;
    MatrixDenseView_new<T> * d_X = nullptr;
    VectorDenseView_new<I> * h_X_colpivots = nullptr;
    VectorDenseView_new<I> * h_L_nnzinsubs = nullptr;
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
    size_t wss_tmp_preprocess_linear = 0;
    size_t wss_tmp_preprocess_overlap = 0;
    size_t wss_tmp_perform_linear = 0;
    size_t wss_tmp_perform_overlap = 0;
    std::unique_ptr<AllocatorArena_new> ator_ws_persistent;
    std::unique_ptr<AllocatorArena_new> ator_ws_tmp_linear;
    std::unique_ptr<AllocatorSinglePointer_new> ator_ws_tmp_overlap;
    MatrixCsxData_new<T,I> d_sub_L_sp;
    MatrixDenseView_new<T> d_sub_L_dn;
    MatrixDenseView_new<T> d_sub_X;
    std::unique_ptr<submatrix_dcsx_dcsx<T,I>> op_d_sub_L_sp;
    math::operations::submatrix_dnx_dnx_view<T> op_sub_L_dn;
    math::operations::submatrix_dnx_dnx_view<T> op_sub_X;
    std::unique_ptr<trsm_dcsx_ddny_ddny<T,I>> op_d_trsm_sp;
    std::unique_ptr<trsm_ddnx_ddny<T>> op_d_trsm_dn;
};



}
}
}

#endif /* SRC_GPU_OPERATIONS_AUXILIARY_GPU_TRSM_TRIRHS_CHUNK_SPITRHS_H */
