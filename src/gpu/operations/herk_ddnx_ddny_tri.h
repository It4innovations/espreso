
#ifndef SRC_GPU_OPERATIONS_HERK_DDNX_DDNY_H
#define SRC_GPU_OPERATIONS_HERK_DDNX_DDNY_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_dnblas.h"
#include "math/math.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class herk_ddnx_ddny_tri
{
    // C = alpha * A * Ah + beta * C
    // C = alpha * Ah * A + beta * C
public:
    struct config
    {
        char strategy = '_'; // sTairs, sQuares
        char partition_algorithm = '_'; // Uniform, Minimal work
        int partition_parameter = 0; // depends on algorithm
    };
    using Treal = utils::remove_complex_t<T>;
public:
    herk_ddnx_ddny_tri() = default;
    herk_ddnx_ddny_tri(const herk_ddnx_ddny_tri &) = delete;
    herk_ddnx_ddny_tri(herk_ddnx_ddny_tri &&) = delete;
    herk_ddnx_ddny_tri & operator=(const herk_ddnx_ddny_tri &) = delete;
    herk_ddnx_ddny_tri & operator=(herk_ddnx_ddny_tri &&) = delete;
    virtual ~herk_ddnx_ddny_tri() = default;
public:
    void set_config(config cfg_);
    void set_handles(gpu::mgm::queue q_, gpu::dnblas::handle dnblas_handle_);
    void set_matrix_d_A(MatrixDenseView_new<T> * d_A_);
    void set_matrix_d_C(MatrixDenseView_new<T> * d_C_);
    void set_h_A_pattern(MatrixCsxView<T,I> * h_A_pattern_);
    void set_coefficients(Treal alpha_, Treal beta_);
    void set_mode(math::herk_mode mode_);
    void setup();
    size_t get_wss_tmp_perform();
    void perform_submit(void * ws_tmp);
private:
    config cfg;
    gpu::mgm::queue q;
    gpu::dnblas::handle handle_dnblas;
    MatrixDenseView_new<T> * d_A = nullptr;
    MatrixDenseView_new<T> * d_C = nullptr;
    MatrixCsxView<T,I> * A_pattern_host = nullptr;
    Treal alpha = Treal{1};
    Treal beta = Treal{0};
    math::herk_mode mode;
    size_t wss_tmp_perform = 0;
    bool called_set_handles = false;
    bool called_set_mode = false;
    bool called_setup = false;
private:
    size_t wss_tmp_peform_linear = 0;
    size_t wss_tmp_peform_overlap = 0;
    std::unique_ptr<AllocatorArena_new> ator_ws_tmp_linear;
    std::unique_ptr<AllocatorSinglePointer_new> ator_ws_tmp_overlap;
    MatrixDenseView_new<T> d_A_reordered;
    MatrixDenseView_new<T> * d_A_left;
    MatrixDenseView_new<T> * d_A_top;
    size_t n;
    size_t k;
    VectorDenseView_new<I> h_A_pivots;
    VectorDenseView_new<I> h_A_trails;
    VectorDenseView_new<size_t> partition;
    size_t num_chunks;
    std::vector<gpu_herk_tri_chunk_stairs<T,I>> op_chunks_stairs;
    std::vector<gpu_herk_tri_chunk_squares<T,I>> op_chunks_squares;
    MatrixDenseView_new<T> d_A_top_dummy;
    herk_ddnx_ddny<T> op_scale_C;
};



}
}
}

#endif /* SRC_GPU_OPERATIONS_HERK_DDNX_DDNY_H */
