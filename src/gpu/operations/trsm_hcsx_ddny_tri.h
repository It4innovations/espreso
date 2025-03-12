
#ifndef SRC_GPU_OPERATIONS_TRSM_HCSX_DDNY_DDNY_TRI_H
#define SRC_GPU_OPERATIONS_TRSM_HCSX_DDNY_DDNY_TRI_H

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
class trsm_hcsx_ddny_tri
{
    // solve A * X = B
public:
    struct config
    {
        char strategy = '_'; // split Factor, split Rhs
    };
public:
    trsm_hcsx_ddny_tri() = default;
    trsm_hcsx_ddny_tri(const trsm_hcsx_ddny_tri &) = delete;
    trsm_hcsx_ddny_tri(trsm_hcsx_ddny_tri &&) = delete;
    trsm_hcsx_ddny_tri & operator=(const trsm_hcsx_ddny_tri &) = delete;
    trsm_hcsx_ddny_tri & operator=(trsm_hcsx_ddny_tri &&) = delete;
    ~trsm_hcsx_ddny_tri() = default;
public:
    void set_handles(gpu::mgm::queue q_, gpu::spblas::handle spblas_handle_, gpu::dnblas::handle dnblas_handle_);
    void set_matrix_L(MatrixCsxView_new<T,I> L);
    void set_matrix_X(MatrixDenseView_new<T> X);
    void calc_X_pattern(MatrixCsxView<T,I> X_pattern_host);
    void setup();
    size_t get_wss_internal();
    size_t get_wss_persistent();
    size_t get_wss_tmp_preprocess();
    size_t get_wss_tmp_perform();
    void set_ws_persistent(void * ws_persistent_);
    void preprocess_submit(void * ws_tmp);
    void update_submit();
    void perform_submit(void * ws_tmp);
private:
    config cfg;
    gpu::mgm::queue q;
    gpu::spblas::handle handle_spblas;
    gpu::dnblas::handle handle_dnblas;
    MatrixCsxView_new<T,I> L;
    MatrixDenseView_new<T> X;
    void * ws_persistent = nullptr;
    size_t wss_internal = 0;
    size_t wss_persistent = 0;
    size_t wss_tmp_preprocess = 0;
    size_t wss_tmp_perform = 0;
    bool called_set_config = false;
    bool called_set_handles = false;
    bool called_set_L = false;
    bool called_set_X = false;
    bool called_calc_X_pattern = false;
    bool called_setup = false;
    bool called_preprocess = false;
private:
    VectorDenseView_new<I> X_colpivots;
    VectorDenseView_new<I> X_rowtrails;
    VectorDenseView_new<I> partition;
};




}
}
}



#endif /* SRC_GPU_OPERATIONS_TRSM_HCSX_DDNY_DDNY_TRI_H */
