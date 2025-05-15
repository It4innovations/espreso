
#ifndef SRC_GPU_OPERATIONS_AUXILIARY_TRSM_HCSX_DDNY_DDNY_TRI_SPLITRHS_H
#define SRC_GPU_OPERATIONS_AUXILIARY_TRSM_HCSX_DDNY_DDNY_TRI_SPLITRHS_H

#include "math/primitives_new/allocator_new.h"
#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/vector_dense_view_new.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_spblas.h"
#include "gpu/gpu_dnblas.h"
#include "gpu/operations/auxiliary/gpu_trsm_trirhs_chunk_splitrhs.h"
#include "math/operations/convert_csx_csy_map.h"
#include "gpu/operations/submatrix_dcsx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class trsm_hcsx_ddny_tri_splitrhs
{
    // solve A * X = B
public:
    struct config
    {
        struct {
            char algorithm = '_'; // Uniform, Minimum work
            int parameter = 0; // meaning depends on algorithm
        } partition;
        char factor_order_sp = '_'; // Rowmajor, Colmajor
        char factor_order_dn = '_'; // Rowmajor, Colmajor
        char spdn_criteria = '_'; // Sparse only, Dense only, fraction of Chunks is sparse, fraction of siZe is sparse, densiTy of factor part
        double spdn_param = 0;
    };
public:
    trsm_hcsx_ddny_tri_splitrhs() = default;
    trsm_hcsx_ddny_tri_splitrhs(const trsm_hcsx_ddny_tri_splitrhs &) = delete;
    trsm_hcsx_ddny_tri_splitrhs(trsm_hcsx_ddny_tri_splitrhs &&) = default;
    trsm_hcsx_ddny_tri_splitrhs & operator=(const trsm_hcsx_ddny_tri_splitrhs &) = delete;
    trsm_hcsx_ddny_tri_splitrhs & operator=(trsm_hcsx_ddny_tri_splitrhs &&) = default;
    ~trsm_hcsx_ddny_tri_splitrhs() = default;
public:
    void set_config(config cfg_);
    void set_handles(gpu::mgm::queue q_, gpu::spblas::handle spblas_handle_, gpu::dnblas::handle dnblas_handle_);
    void set_matrix_h_L(MatrixCsxView_new<T,I> * h_L_);
    void set_matrix_d_X(MatrixDenseView_new<T> * d_X_);
    void set_h_X_pattern(MatrixCsxView_new<T,I> * h_X_pattern_);
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
    gpu::mgm::queue q;
    gpu::spblas::handle handle_spblas;
    gpu::dnblas::handle handle_dnblas;
    MatrixCsxView_new<T,I> * h_L = nullptr;
    MatrixDenseView_new<T> * d_X = nullptr;
    MatrixCsxView_new<T,I> * h_X_pattern = nullptr;
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
    std::unique_ptr<AllocatorArena_new> ator_ws_persistent;
    std::unique_ptr<AllocatorArena_new> ator_ws_tmp_linear;
    std::unique_ptr<AllocatorSinglePointer_new> ator_ws_tmp_overlap;
    size_t wss_tmp_preprocess_linear = 0;
    size_t wss_tmp_preprocess_overlap = 0;
    size_t wss_tmp_perform_linear = 0;
    size_t wss_tmp_perform_overlap = 0;
    VectorDenseData_new<I> h_X_colpivots;
    VectorDenseData_new<I> h_X_rowtrails;
    VectorDenseData_new<size_t> partition;
    VectorDenseData_new<I> h_L_nnzinsubs_splitrhs;
    size_t num_chunks = 0;
    size_t first_dense_chunk;
    std::vector<gpu_trsm_trirhs_chunk_splitrhs<T,I>> ops_chunks;
    MatrixCsxData_new<T,I> h_L_reordered;
    MatrixCsxView_new<T,I> * h_L_to_use = nullptr;
    MatrixCsxData_new<T,I> d_L_sp;
    MatrixDenseData_new<T> d_L_dn;
    math::operations::convert_csx_csy_map<T,I> op_h_L_reorder;
    std::unique_ptr<submatrix_dcsx_ddny<T,I>> op_d_sub_L_sp2dn;
};




}
}
}



#endif /* SRC_GPU_OPERATIONS_AUXILIARY_TRSM_HCSX_DDNY_DDNY_TRI_SPLITRHS_H */
