
#include "math/operations/herk_dnx_dny_tri.h"

#include "math/operations/pivots_trails_csx.h"
#include "math/operations/auxiliary/tri_partition_herk.h"
#include "math/operations/herk_dnx_dny.h"
#include "math/operations/gemm_dnx_dny_dnz.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
herk_dnx_dny_tri<T,I>::~herk_dnx_dny_tri()
{
    finalize();
}



template<typename T, typename I>
void herk_dnx_dny_tri<T,I>::set_config(config cfg_)
{
    cfg = cfg_;

    config_set = true;
}



template<typename T, typename I>
void herk_dnx_dny_tri<T,I>::set_matrix_A(MatrixDenseView_new<T> * A_)
{
    A = A_;
}



template<typename T, typename I>
void herk_dnx_dny_tri<T,I>::set_matrix_C(MatrixDenseView_new<T> * C_)
{
    C = C_;
}



template<typename T, typename I>
void herk_dnx_dny_tri<T,I>::set_coefficients(Treal alpha_, Treal beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void herk_dnx_dny_tri<T,I>::set_mode(blas::herk_mode mode_)
{
    mode = mode_;

    mode_set = true;
}



template<typename T, typename I>
void herk_dnx_dny_tri<T,I>::calc_A_pattern(MatrixCsxView_new<T,I> & A_pattern)
{
    stacktimer::push("herk_dnx_dny_tri::calc_A_pattern");

    if(!mode_set) eslog::error("mode is not set\n");
    if(pattern_set) eslog::error("A pattern was already set\n");

    if(mode == blas::herk_mode::AhA) {
        A_pivots.set(A->ncols, AllocatorCPU_new::get_singleton());
        A_pivots.alloc();
        pivots_trails_csx<T,I>::do_all(&A_pattern, &A_pivots, 'C', 'P', 'B');

        A_trails.set(A->nrows, AllocatorCPU_new::get_singleton());
        A_trails.alloc();
        pivots_trails_csx<T,I>::do_all(&A_pattern, &A_trails, 'R', 'T', 'F');
    }
    if(mode == blas::herk_mode::AAh) {
        A_pivots.set(A->nrows, AllocatorCPU_new::get_singleton());
        A_pivots.alloc();
        pivots_trails_csx<T,I>::do_all(&A_pattern, &A_pivots, 'R', 'P', 'B');

        A_trails.set(A->ncols, AllocatorCPU_new::get_singleton());
        A_trails.alloc();
        pivots_trails_csx<T,I>::do_all(&A_pattern, &A_trails, 'C', 'T', 'F');
    }

    stacktimer::pop();

    pattern_set = true;
}



template<typename T, typename I>
void herk_dnx_dny_tri<T,I>::preprocess()
{
    stacktimer::push("herk_dnx_dny_tri::preprocess");

    if(!config_set) eslog::error("config is not set\n");
    if(!mode_set) eslog::error("mode is not set\n");
    if(!pattern_set) eslog::error("A pattern is not set\n");

    for(size_t i = 1; i < A_pivots.size; i++) {
        if(A_pivots.vals[i-1] > A_pivots.vals[i]) {
            eslog::error("A does not have proper triangular structure\n");
        }
    }
    for(size_t i = 1; i < A_trails.size; i++) {
        if(A_trails.vals[i-1] > A_trails.vals[i]) {
            eslog::error("A does not have proper triangular structure\n");
        }
    }

    n = A_pivots.size;
    k = A_trails.size;

    tri_partition_herk partitioner;
    char partition_direction = '_';
    if(cfg.strategy == 'T') partition_direction = 'N';
    if(cfg.strategy == 'Q') partition_direction = 'K';
    partitioner.set_config(cfg.partition_algorithm, partition_direction, cfg.partition_parameter, cfg.strategy);
    partitioner.set_system(n, k);
    partitioner.set_output_partition(&partition);
    partitioner.setup();
    num_chunks = partitioner.get_num_chunks();
    partition.set(num_chunks + 1, AllocatorCPU_new::get_singleton());
    partition.alloc();
    partitioner.perform();

    stacktimer::pop();

    preproces_called = true;
}



template<typename T, typename I>
void herk_dnx_dny_tri<T,I>::perform()
{
    stacktimer::push("herk_dnx_dny_tri::perform");

    if(!preproces_called) eslog::error("preprocess was not called\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(C->nrows != C->ncols) eslog::error("C must be square\n");
    if(C->prop.uplo != 'U' && C->prop.uplo != 'L') eslog::error("wrong C uplo\n");
    if(A->ncols != C->ncols) eslog::error("incompatible matrices\n");

    if(mode == blas::herk_mode::AAh) {
        perform_AAh();
    }
    if(mode == blas::herk_mode::AhA) {
        perform_AhA();
    }

    stacktimer::pop();
}



template<typename T, typename I>
void herk_dnx_dny_tri<T,I>::perform_AhA()
{
    MatrixDenseView_new<T> sub_A_trivial = A->get_submatrix_view(0, 0, 0, n);
    herk_dnx_dny<T>::do_all(&sub_A_trivial, C, blas::herk_mode::AhA, Treal{0}, beta);

    for(size_t ch = 0; ch < num_chunks; ch++) {
        if(cfg.strategy == 'T') {
            size_t n_start = partition.vals[ch];
            size_t n_end = partition.vals[ch+1];
            size_t k_start = A_pivots.vals[n_start];

            if(C->prop.uplo == 'L' && ch > 0) {
                MatrixDenseView_new<T> sub_A_left = A->get_transposed_reordered_view().get_submatrix_view(n_start, n_end, k_start, k);
                MatrixDenseView_new<T> sub_A_top = A->get_submatrix_view(k_start, k, 0, n_start);
                MatrixDenseView_new<T> sub_C = C->get_submatrix_view(n_start, n_end, 0, n_start);
                sub_C.prop.uplo = C->prop.uplo;
                gemm_dnx_dny_dnz<T>::do_all(&sub_A_left, &sub_A_top, &sub_C, alpha, T{1}, true, false);
            }
            if(C->prop.uplo == 'U' && ch > 0) {
                MatrixDenseView_new<T> sub_A_left = A->get_transposed_reordered_view().get_submatrix_view(0, n_start, k_start, k);
                MatrixDenseView_new<T> sub_A_top = A->get_submatrix_view(k_start, k, n_start, n_end);
                MatrixDenseView_new<T> sub_C = C->get_submatrix_view(0, n_start, n_start, n_end);
                sub_C.prop.uplo = C->prop.uplo;
                gemm_dnx_dny_dnz<T>::do_all(&sub_A_left, &sub_A_top, &sub_C, alpha, T{1}, true, false);
            }
            {
                MatrixDenseView_new<T> sub_A = A->get_submatrix_view(k_start, k, n_start, n_end);
                MatrixDenseView_new<T> sub_C = C->get_submatrix_view(n_start, n_end, n_start, n_end);
                sub_C.prop.uplo = C->prop.uplo;
                herk_dnx_dny<T>::do_all(&sub_A, &sub_C, mode, alpha, Treal{1});
            }
        }
        if(cfg.strategy == 'Q') {
            size_t k_start = partition.vals[ch];
            size_t k_end = partition.vals[ch+1];
            size_t n_end = A_trails.vals[k_end - 1] + 1;

            MatrixDenseView_new<T> sub_A = A->get_submatrix_view(k_start, k_end, 0, n_end);
            MatrixDenseView_new<T> sub_C = C->get_submatrix_view(0, n_end, 0, n_end);
            sub_C.prop.uplo = C->prop.uplo;
            herk_dnx_dny<T>::do_all(&sub_A, &sub_C, mode, alpha, Treal{1});
        }
    }
}



template<typename T, typename I>
void herk_dnx_dny_tri<T,I>::perform_AAh()
{
    MatrixDenseView_new<T> sub_A_trivial = A->get_submatrix_view(0, n, 0, 0);
    herk_dnx_dny<T>::do_all(&sub_A_trivial, C, blas::herk_mode::AAh, Treal{0}, beta);

    for(size_t ch = 0; ch < num_chunks; ch++) {
        if(cfg.strategy == 'T') {
            size_t n_start = partition.vals[ch];
            size_t n_end = partition.vals[ch+1];
            size_t k_start = A_pivots.vals[n_start];

            if(C->prop.uplo == 'L' && ch > 0) {
                MatrixDenseView_new<T> sub_A_left = A->get_submatrix_view(n_start, n_end, k_start, k);
                MatrixDenseView_new<T> sub_A_top = A->get_transposed_reordered_view().get_submatrix_view(k_start, k, 0, n_start);
                MatrixDenseView_new<T> sub_C = C->get_submatrix_view(n_start, n_end, 0, n_start);
                sub_C.prop.uplo = C->prop.uplo;
                gemm_dnx_dny_dnz<T>::do_all(&sub_A_left, &sub_A_top, &sub_C, alpha, T{1}, false, true);
            }
            if(C->prop.uplo == 'U' && ch > 0) {
                MatrixDenseView_new<T> sub_A_left = A->get_submatrix_view(0, n_start, k_start, k);
                MatrixDenseView_new<T> sub_A_top = A->get_transposed_reordered_view().get_submatrix_view(k_start, k, n_start, n_end);
                MatrixDenseView_new<T> sub_C = C->get_submatrix_view(0, n_start, n_start, n_end);
                sub_C.prop.uplo = C->prop.uplo;
                gemm_dnx_dny_dnz<T>::do_all(&sub_A_left, &sub_A_top, &sub_C, alpha, T{1}, false, true);
            }
            {
                MatrixDenseView_new<T> sub_A = A->get_submatrix_view(n_start, n_end, k_start, k);
                MatrixDenseView_new<T> sub_C = C->get_submatrix_view(n_start, n_end, n_start, n_end);
                sub_C.prop.uplo = C->prop.uplo;
                herk_dnx_dny<T>::do_all(&sub_A, &sub_C, mode, alpha, Treal{1});
            }
        }
        if(cfg.strategy == 'Q') {
            size_t k_start = partition.vals[ch];
            size_t k_end = partition.vals[ch+1];
            size_t n_end = A_trails.vals[k_end - 1] + 1;

            MatrixDenseView_new<T> sub_A = A->get_submatrix_view(0, n_end, k_start, k_end);
            MatrixDenseView_new<T> sub_C = C->get_submatrix_view(0, n_end, 0, n_end);
            sub_C.prop.uplo = C->prop.uplo;
            herk_dnx_dny<T>::do_all(&sub_A, &sub_C, mode, alpha, Treal{1});
        }
    }
}



template<typename T, typename I>
void herk_dnx_dny_tri<T,I>::finalize()
{
    if(preproces_called) {
        partition.clear();
    }
    preproces_called = false;

    if(pattern_set) {
        A_pivots.clear();
        A_trails.clear();
    }
    pattern_set = false;
}



#define INSTANTIATE_T_I(T,I) \
template class herk_dnx_dny_tri<T,I>; \

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T,int32_t) \
    /* INSTANTIATE_T_I(T,int64_t) */

        #define INSTANTIATE \
        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        INSTANTIATE_T(std::complex<double>)

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}
