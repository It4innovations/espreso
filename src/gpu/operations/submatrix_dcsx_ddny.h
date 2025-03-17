
#ifndef SRC_GPU_OPERATIONS_SUBMATRIX_DCSX_DDNY_H
#define SRC_GPU_OPERATIONS_SUBMATRIX_DCSX_DDNY_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"
#include "gpu/gpu_management.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class submatrix_dcsx_ddny
{
    // does not respect uplo and diag
protected:
    submatrix_dcsx_ddny() = default;
public:
    submatrix_dcsx_ddny(const submatrix_dcsx_ddny &) = delete;
    submatrix_dcsx_ddny(submatrix_dcsx_ddny &&) = delete;
    submatrix_dcsx_ddny & operator=(const submatrix_dcsx_ddny &) = delete;
    submatrix_dcsx_ddny & operator=(submatrix_dcsx_ddny &&) = delete;
    virtual ~submatrix_dcsx_ddny() = default;
public:
    static std::unique_ptr<submatrix_dcsx_ddny<T,I>> make();
public:
    void set_handles(gpu::mgm::queue q_);
    void set_bounds(size_t row_start_, size_t row_end_, size_t col_start_, size_t col_end_);
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * M_dst_);
    void setup();
    size_t get_wss_internal();
    size_t get_wss_persistent();
    size_t get_wss_tmp_preprocess();
    size_t get_wss_tmp_perform();
    void set_ws_persistent(void * ws_persistent_);
    void preprocess_submit(void * ws_tmp);
    void perform_submit(void * ws_tmp);
protected:
    gpu::mgm::queue q;
    size_t row_start = 0;
    size_t row_end = 0;
    size_t col_start = 0;
    size_t col_end = 0;
    size_t primary_start = 0;
    size_t primary_end = 0;
    size_t secdary_start = 0;
    size_t secdary_end = 0;
    MatrixCsxView_new<T,I> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
    void * ws_persistent = nullptr;
    size_t wss_internal = 0;
    size_t wss_persistent = 0;
    size_t wss_tmp_preprocess = 0;
    size_t wss_tmp_perform = 0;
    bool called_set_handles = false;
    bool called_set_bounds = false;
    bool called_setup = false;
    bool called_preprocess = false;
protected:
    virtual void internal_setup() {}
    virtual void internal_preprocess(void * ws_tmp) {}
    virtual void internal_perform(void * ws_tmp) {}
};



}
}
}



#endif /* SRC_GPU_OPERATIONS_SUBMATRIX_DCSX_DDNY_H */
