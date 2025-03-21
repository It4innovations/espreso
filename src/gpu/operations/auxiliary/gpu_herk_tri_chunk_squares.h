
#ifndef SRC_GPU_OPERATIONS_AUXILIARY_GPU_HERK_TRI_CHUNK_SQUARES_H
#define SRC_GPU_OPERATIONS_AUXILIARY_GPU_HERK_TRI_CHUNK_SQUARES_H

#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/operations/submatrix_dnx_dnx_view.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_dnblas.h"
#include "gpu/operations/herk_ddnx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class gpu_herk_tri_chunk_squares
{
public:
    using Treal = utils::remove_complex_t<T>;
public:
    gpu_herk_tri_chunk_squares() = default;
    gpu_herk_tri_chunk_squares(const gpu_herk_tri_chunk_squares &) = delete;
    gpu_herk_tri_chunk_squares(gpu_herk_tri_chunk_squares &&) = default;
    gpu_herk_tri_chunk_squares & operator=(const gpu_herk_tri_chunk_squares &) = delete;
    gpu_herk_tri_chunk_squares & operator=(gpu_herk_tri_chunk_squares &&) = default;
    ~gpu_herk_tri_chunk_squares() = default;
public:
    void set_range(size_t k_start_, size_t k_end_);
    void set_handles(gpu::mgm::queue q_, gpu::dnblas::handle handle_dnblas_);
    void set_matrix_d_A_left(MatrixDenseView_new<T> * d_A_left_);
    void set_matrix_d_A_top(MatrixDenseView_new<T> * d_A_top_);
    void set_matrix_d_C(MatrixDenseView_new<T> * d_C_);
    void set_h_A_trails(VectorDenseView_new<I> * h_A_trails_);
    void set_coefficients(Treal alpha_);
    void setup();
    size_t get_wss_tmp_perform();
    void perform_submit(void * ws_tmp);
private:
    size_t k_start = 0;
    size_t k_end = 0;
    gpu::mgm::queue q;
    gpu::dnblas::handle handle_dnblas;
    MatrixDenseView_new<T> * d_A_left = nullptr;
    MatrixDenseView_new<T> * d_A_top = nullptr;
    MatrixDenseView_new<T> * d_C = nullptr;
    VectorDenseView_new<I> * h_A_trails = nullptr;
    Treal alpha = Treal{1};
    size_t wss_tmp_perform = 0;
    bool called_set_range = false;
    bool called_set_handles = false;
    bool called_setup = false;
private:
    size_t k_size = 0;
    size_t n_end = 0;
    size_t n_size = 0;
    MatrixDenseView_new<T> d_sub_C;
    MatrixDenseView_new<T> d_sub_A_top;
    math::operations::submatrix_dnx_dnx_view<T> op_sub_C;
    math::operations::submatrix_dnx_dnx_view<T> op_sub_A_top;
    std::unique_ptr<herk_ddnx_ddny<T>> op_herk;
};



}
}
}

#endif /* SRC_GPU_OPERATIONS_AUXILIARY_GPU_HERK_TRI_CHUNK_SQUARES_H */
