
#ifndef SRC_GPU_OPERATIONS_SUPERMATRIX_DDNX_DDNX_NONCONTIG_H
#define SRC_GPU_OPERATIONS_SUPERMATRIX_DDNX_DDNX_NONCONTIG_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/vector_dense_view_new.h"
#include "gpu/gpu_management.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class supermatrix_ddnx_ddnx_noncontig
{
    // ignores uplo
public:
    enum struct mode {
        assign,
        accumulate,
        accumulate_atomic,
    };
protected:
    supermatrix_ddnx_ddnx_noncontig() = default;
public:
    supermatrix_ddnx_ddnx_noncontig(const supermatrix_ddnx_ddnx_noncontig &) = delete;
    supermatrix_ddnx_ddnx_noncontig(supermatrix_ddnx_ddnx_noncontig &&) = delete;
    supermatrix_ddnx_ddnx_noncontig & operator=(const supermatrix_ddnx_ddnx_noncontig &) = delete;
    supermatrix_ddnx_ddnx_noncontig & operator=(supermatrix_ddnx_ddnx_noncontig &&) = delete;
    virtual ~supermatrix_ddnx_ddnx_noncontig() = default;
public:
    static std::unique_ptr<supermatrix_ddnx_ddnx_noncontig<T,I>> make();
public:
    void set_handles(gpu::mgm::queue q_);
    void set_matrix_src(MatrixDenseView_new<T> * d_M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * d_M_dst_);
    void set_row_map(VectorDenseView_new<I> * d_row_map_);
    void set_col_map(VectorDenseView_new<I> * d_col_map_);
    void set_mode(mode mode_val_);
    void perform_submit();
    static void submit_all(gpu::mgm::queue q, MatrixDenseView_new<T> * d_M_src, MatrixDenseView_new<T> * d_M_dst, VectorDenseView_new<I> * d_row_map, VectorDenseView_new<I> * d_col_map, mode mode_val = mode::assign);
protected:
    gpu::mgm::queue q;
    MatrixDenseView_new<T> * d_M_src = nullptr;
    MatrixDenseView_new<T> * d_M_dst = nullptr;
    VectorDenseView_new<I> * d_row_map = nullptr;
    VectorDenseView_new<I> * d_col_map = nullptr;
    mode mode_val = mode::assign;
    bool called_set_handles = false;
    bool called_set_row_map = false;
    bool called_set_col_map = false;
protected:
    virtual void internal_perform() {}
};



}
}
}

#endif /* SRC_GPU_OPERATIONS_SUPERMATRIX_DDNX_DDNX_NONCONTIG_H */
