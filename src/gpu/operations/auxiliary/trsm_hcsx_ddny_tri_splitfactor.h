
#ifndef SRC_GPU_OPERATIONS_AUXILIARY_TRSM_HCSX_DDNY_DDNY_TRI_SPLITFACTOR_H
#define SRC_GPU_OPERATIONS_AUXILIARY_TRSM_HCSX_DDNY_DDNY_TRI_SPLITFACTOR_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/vector_dense_view_new.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_spblas.h"
#include "gpu/gpu_dnblas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class trsm_hcsx_ddny_tri_splitfactor
{
    // solve A * X = B
public:
    struct config
    {
        struct {
            char algorithm = '_'; // Uniform, Minimum work
            int parameter = 0; // meaning depends on algorithm
        } partition;
        char trsm_factor_spdn = '_'; // Sparse, Dense
        char trsm_factor_order = '_'; // Rowmajor, Colmajor
        char gemm_factor_prune = '_'; // No, Rows only, Cols only, All
        char gemm_factor_order_sp = '_'; // Rowmajor, Colmajor
        char gemm_factor_order_dn = '_'; // Rowmajor, Colmajor
        char gemm_spdn_criteria = '_'; // Sparse only, Dense only, fraction of Chunks is sparse, fraction of siZe is sparse, densiTy of factor part
        double gemm_spdn_param = 0;
    };
public:
    trsm_hcsx_ddny_tri_splitfactor() = default;
    trsm_hcsx_ddny_tri_splitfactor(const trsm_hcsx_ddny_tri_splitfactor &) = delete;
    trsm_hcsx_ddny_tri_splitfactor(trsm_hcsx_ddny_tri_splitfactor &&) = delete;
    trsm_hcsx_ddny_tri_splitfactor & operator=(const trsm_hcsx_ddny_tri_splitfactor &) = delete;
    trsm_hcsx_ddny_tri_splitfactor & operator=(trsm_hcsx_ddny_tri_splitfactor &&) = delete;
    ~trsm_hcsx_ddny_tri_splitfactor() = default;
public:
    void set_config(config cfg_);
    void set_handles(gpu::mgm::queue q_, gpu::spblas::handle spblas_handle_, gpu::dnblas::handle dnblas_handle_);
    void set_matrix_h_L(MatrixCsxView_new<T,I> * h_L_);
    void set_matrix_d_X(MatrixDenseView_new<T> * d_X_);
    void calc_X_pattern(MatrixCsxView<T,I> & X_pattern_host);
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
    void * ws_persistent = nullptr;
    size_t wss_internal = 0;
    size_t wss_persistent = 0;
    size_t wss_tmp_preprocess = 0;
    size_t wss_tmp_perform = 0;
    bool called_set_config = false;
    bool called_set_handles = false;
    bool called_calc_X_pattern = false;
    bool called_setup = false;
    bool called_preprocess = false;
private:
    size_t wss_tmp_preprocess_linear = 0;
    size_t wss_tmp_preprocess_overlap = 0;
    size_t wss_tmp_peform_linear = 0;
    size_t wss_tmp_peform_overlap = 0;
    std::unique_ptr<AllocatorArena_new> ator_ws_persistent;
    std::unique_ptr<AllocatorArena_new> ator_ws_tmp_linear;
    std::unique_ptr<AllocatorSinglePointer_new> ator_ws_tmp_overlap;
    VectorDenseView_new<I> h_X_colpivots;
    VectorDenseView_new<I> h_X_rowtrails;
    VectorDenseView_new<I> h_partition;
    size_t num_chunks = 0;
    std::vector<gpu_trsm_trirhs_chunk_splitfactor<T,I>> ops_chunks;
};




}
}
}



#endif /* SRC_GPU_OPERATIONS_AUXILIARY_TRSM_HCSX_DDNY_DDNY_TRI_SPLITFACTOR_H */
