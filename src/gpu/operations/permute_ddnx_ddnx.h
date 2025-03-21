
#ifndef SRC_GPU_OPERATIONS_PERMUTE_DDNX_DDNX_H
#define SRC_GPU_OPERATIONS_PERMUTE_DDNX_DDNX_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/permutation_view_new.h"
#include "gpu/gpu_management.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class permute_ddnx_ddnx
{
    // does not respect uplo
protected:
    permute_ddnx_ddnx() = default;
public:
    permute_ddnx_ddnx(const permute_ddnx_ddnx &) = delete;
    permute_ddnx_ddnx(permute_ddnx_ddnx &&) = delete;
    permute_ddnx_ddnx & operator=(const permute_ddnx_ddnx &) = delete;
    permute_ddnx_ddnx & operator=(permute_ddnx_ddnx &&) = delete;
    virtual ~permute_ddnx_ddnx() = default;
public:
    static std::unique_ptr<permute_ddnx_ddnx<T,I>> make();
public:
    void set_handles(gpu::mgm::queue q_);
    void set_matrix_src(MatrixDenseView_new<T> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * M_dst_);
    void set_perm_rows(PermutationView_new<I> * perm_rows_);
    void set_perm_cols(PermutationView_new<I> * perm_cols_);
    void setup();
    size_t get_wss_tmp_perform();
    void perform_submit(void * ws_tmp);
protected:
    gpu::mgm::queue q;
    MatrixDenseView_new<T> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
    PermutationView_new<I> * perm_rows = nullptr;
    PermutationView_new<I> * perm_cols = nullptr;
    size_t wss_tmp_perform = 0;
    bool called_set_handles = false;
    bool called_setup = false;
    PermutationView_new<I> * perm_primary = nullptr;
    PermutationView_new<I> * perm_secdary = nullptr;
protected:
    virtual void internal_setup() {}
    virtual void internal_perform(void * /*ws_tmp*/) {}
};



}
}
}

#endif /* SRC_GPU_OPERATIONS_PERMUTE_DDNX_DDNX_H */
