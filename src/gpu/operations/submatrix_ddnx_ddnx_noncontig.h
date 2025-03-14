
#ifndef SRC_GPU_OPERATIONS_SUBMATRIX_DDNX_DDNX_NONCONTIG_H
#define SRC_GPU_OPERATIONS_SUBMATRIX_DDNX_DDNX_NONCONTIG_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/vector_dense_view_new.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class submatrix_ddnx_ddnx_noncontig
{
protected:
    submatrix_ddnx_ddnx_noncontig() = default;
public:
    submatrix_ddnx_ddnx_noncontig(const submatrix_ddnx_ddnx_noncontig &) = delete;
    submatrix_ddnx_ddnx_noncontig(submatrix_ddnx_ddnx_noncontig &&) = delete;
    submatrix_ddnx_ddnx_noncontig & operator=(const submatrix_ddnx_ddnx_noncontig &) = delete;
    submatrix_ddnx_ddnx_noncontig & operator=(submatrix_ddnx_ddnx_noncontig &&) = delete;
    ~submatrix_ddnx_ddnx_noncontig() = default;
public:
    static std::unique_ptr<submatrix_ddnx_ddnx_noncontig<T,I>> make();
public:
    void set_handles(gpu::mgm::queue q_);
    void set_matrix_src(MatrixDenseView_new<T> * d_M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * d_M_dst_);
    void set_row_map(VectorDenseView_new<I> * d_row_map_);
    void set_col_map(VectorDenseView_new<I> * d_col_map_);
    void perform_submit();
    static void submit_all(gpu::mgm::queue q, MatrixDenseView_new<T> * d_M_src, MatrixDenseView_new<T> * d_M_dst, VectorDenseView_new<I> * d_row_map, VectorDenseView_new<I> * d_col_map);
protected:
    gpu::mgm::queue q;
    MatrixDenseView_new<T> * d_M_src = nullptr;
    MatrixDenseView_new<T> * d_M_dst = nullptr;
    VectorDenseView_new<I> * d_row_map = nullptr;
    VectorDenseView_new<I> * d_col_map = nullptr;
    bool called_set_handles = false;
    bool called_set_row_map = false;
    bool called_set_col_map = false;
protected:
    virtual void internal_perform();
};



}
}
}

#endif /* SRC_GPU_OPERATIONS_SUBMATRIX_DDNX_DDNX_NONCONTIG_H */
