
#include "math/operations/herk_dnx_dny_tri.h"



template<typename T, typename I>
herk_dnx_dny_tri<T>::~herk_dnx_dny_tri()
{
    finalize();
}



template<typename T, typename I>
void herk_dnx_dny_tri<T>::set_config(config cfg_)
{
    cfg = cfg_;

    config_set = true;
}



template<typename T, typename I>
void herk_dnx_dny_tri<T>::set_matrix_A(MatrixDenseView_new<T> * A_)
{
    A = A_;
}



template<typename T, typename I>
void herk_dnx_dny_tri<T>::set_matrix_C(MatrixDenseView_new<T> * C_)
{
    C = C_;
}



template<typename T, typename I>
void herk_dnx_dny_tri<T>::set_coefficients(T alpha_, T beta_)
{
    alpha = alpha_;
    beta = beta_;
}



template<typename T, typename I>
void herk_dnx_dny_tri<T>::set_mode(blas::herk_mode mode_)
{
    mode = mode_;

    mode_set = true;
}



template<typename T, typename I>
void herk_dnx_dny_tri<T>::set_A_pattern(MatrixCsxView_new<T,I> & A_pattern)
{
    if(!mode_set) eslog::error("mode is not set\n");
    if(pattern_set) eslog::error("A pattern was already set\n");

    if(mode == blas::herk_mode::AtA) {
        A_pivots.set(A.ncols, AllocatorCPU_new::get_singleton());
        A_pivots.alloc();
        pivots_trails_csx<T,I>::do_all(&A_pattern, &A_pivots, 'C', 'P');

        A_trails.set(A.nrows, AllocatorCPU_new::get_singleton());
        A_trails.alloc();
        pivots_trails_csx<T,I>::do_all(&A_pattern, &A_trails, 'R', 'T');
    }
    if(mode == blas::herk_mode::AAt) {
        A_pivots.set(A.nrows, AllocatorCPU_new::get_singleton());
        A_pivots.alloc();
        pivots_trails_csx<T,I>::do_all(&A_pattern, &A_pivots, 'R', 'P');

        A_trails.set(A.ncols, AllocatorCPU_new::get_singleton());
        A_trails.alloc();
        pivots_trails_csx<T,I>::do_all(&A_pattern, &A_trails, 'C', 'T');
    }

    pattern_set = true;
}



template<typename T, typename I>
void herk_dnx_dny_tri<T>::preprocess()
{
    if(!config_set) eslog::error("config is not set\n");
    if(!mode_set) eslog::error("mode is not set\n");
    if(!pattern_set) eslog::error("A pattern is not set\n");
    if(A == nullptr) eslog::error("matrix A is not set\n");
    if(C == nullptr) eslog::error("matrix C is not set\n");
    if(A_colpivots == nullptr) eslog::error("A colpivots are not set\n");
    if(A_rowtrails == nullptr) eslog::error("A rowtrails are not set\n");
    if(C->nrows != C->ncols) eslog::error("C must be square\n");
    if(C->uplo != 'U' && C->uplo != 'L') eslog::error("wrong C uplo\n");
    if(A->ncols != C->ncols) eslog::error("incompatible matrices\n");

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
    if(strategy == 'T') partition_direction = 'N';
    if(strategy == 'Q') partition_direction = 'K';
    partitioner.set_config(cfg.partition.algorithm, partition_direction, cfg.partition.parameter, cfg.strategy);
    partitioner.set_system(n, k);
    partitioner.set_output_partition(partition);
    partitioner.setup();
    num_chunks = partitioner.get_num_chunks();
    partition.set(num_chunks, AllocatorCPU_new::get_singleton());
    partition.alloc();
    partitioner.perform();

    preproces_called = true;
}



template<typename T, typename I>
void herk_dnx_dny_tri<T>::perform()
{
    if(!preproces_called) eslog::error("preprocess was not called\n");

    if(mode == blas::herk_mode::AAh) {
        perform_AAh();
    }
    if(mode == blas::herk_mode::AAh) {
        perform_AAh();
    }    
}



template<typename T, typename I>
void herk_dnx_dny_tri<T>::perform_AhA()
{
    MatrixDenseView_new<T> sub_A_trivial = A->make_submatrix(0, 0, 0, n);
    herk_dnx_dny<T>::do_all(sub_A_trivial, C, blas::herk_mode::AhA, T{0}, beta);

    for(size_t ch = 0; ch < num_chunks; ch++) {
        if(cfg.strategy == 'T') {
            size_t n_start = partition[ch];
            size_t n_end = partition[ch+1];
            size_t k_start = B_pivots[n_start];

            if(C->uplo == 'L' && ch > 0) {
                MatrixDenseView_new<T> sub_A_left = A->make_transposed_reordered().make_submatrix(n_start, n_end, k_start, k);
                MatrixDenseView_new<T> sub_A_top = A->make_submatrix(k_start, k, 0, n_start);
                MatrixDenseView_new<T> sub_C = C->make_submatrix(n_start, n_end, 0, n_start);
                sub_C.prop.uplo = C->prop.uplo;
                gemm_dnx_dny_dnz<T>::do_all(&sub_A_left, &sub_A_top, &sub_C, alpha, T{1}, true, false);
            }
            if(C->uplo == 'U' && ch > 0) {
                MatrixDenseView_new<T> sub_A_left = A->make_transposed_reordered().make_submatrix(0, n_start, k_start, k);
                MatrixDenseView_new<T> sub_A_top = A->make_submatrix(k_start, k, n_start, n_end);
                MatrixDenseView_new<T> sub_C = C->make_submatrix(0, n_start, n_start, n_end);
                sub_C.prop.uplo = C->prop.uplo;
                gemm_dnx_dny_dnz<T>::do_all(&sub_A_left, &sub_A_top, &sub_C, alpha, T{1}, true, false);
            }
            {
                MatrixDenseView_new<T> sub_A = A->make_submatrix(k_start, k, n_start, n_end);
                MatrixDenseView_new<T> sub_C = C->make_submatrix(n_start, n_end, n_start, n_end);
                sub_C.prop.uplo = C->prop.uplo;
                herk_dnx_dny<T>::do_all(&sub_A, &sub_C, mode, alpha, T{1});
            }
        }
        if(cfg.strategy == 'Q') {
            size_t k_start = partition[ch];
            size_t k_end = partition[ch+1];
            size_t n_end = A_trails[k_end - 1] + 1;

            MatrixDenseView_new<T> sub_A = A->make_submatrix(k_start, k_end, n_start, n_end);
            MatrixDenseView_new<T> sub_C = C->make_submatrix(n_start, n_end, n_start, n_end);
            sub_C.prop.uplo = C->prop.uplo;
            herk_dnx_dny<T>::do_all(sub_A, sub_C, mode, alpha, T{1});
        }
    }
}



template<typename T, typename I>
void herk_dnx_dny_tri<T>::perform_AAh()
{
    MatrixDenseView_new<T> sub_A_trivial = A->make_submatrix(0, n, 0, 0);
    herk_dnx_dny<T>::do_all(sub_A_trivial, C, blas::herk_mode::AAh, T{0}, beta);

    for(size_t ch = 0; ch < num_chunks; ch++) {
        if(cfg.strategy == 'T') {
            size_t n_start = partition[ch];
            size_t n_end = partition[ch+1];
            size_t k_start = B_pivots[n_start];

            if(C->uplo == 'L' && ch > 0) {
                MatrixDenseView_new<T> sub_A_left = A->make_submatrix(n_start, n_end, k_start, k);
                MatrixDenseView_new<T> sub_A_top = A->make_transposed_reordered().make_submatrix(k_start, k, 0, n_start);
                MatrixDenseView_new<T> sub_C = C->make_submatrix(n_start, n_end, 0, n_start);
                sub_C.prop.uplo = C->prop.uplo;
                gemm_dnx_dny_dnz<T>::do_all(&sub_A_left, &sub_A_top, &sub_C, alpha, T{1}, false, true);
            }
            if(C->uplo == 'U' && ch > 0) {
                MatrixDenseView_new<T> sub_A_left = A->make_submatrix(0, n_start, k_start, k);
                MatrixDenseView_new<T> sub_A_top = A->make_transposed_reordered().make_submatrix(k_start, k, n_start, n_end);
                MatrixDenseView_new<T> sub_C = C->make_submatrix(0, n_start, n_start, n_end);
                sub_C.prop.uplo = C->prop.uplo;
                gemm_dnx_dny_dnz<T>::do_all(&sub_A_left, &sub_A_top, &sub_C, alpha, T{1}, false, true);
            }
            {
                MatrixDenseView_new<T> sub_A = A->make_submatrix(n_start, n_end, k_start, k);
                MatrixDenseView_new<T> sub_C = C->make_submatrix(n_start, n_end, n_start, n_end);
                sub_C.prop.uplo = C->prop.uplo;
                herk_dnx_dny<T>::do_all(&sub_A, &sub_C, mode, alpha, T{1});
            }
        }
        if(cfg.strategy == 'Q') {
            size_t k_start = partition[ch];
            size_t k_end = partition[ch+1];
            size_t n_end = A_trails[k_end - 1] + 1;

            MatrixDenseView_new<T> sub_A = A->make_submatrix(n_start, n_end, k_start, k_end);
            MatrixDenseView_new<T> sub_C = C->make_submatrix(n_start, n_end, n_start, n_end);
            sub_C.prop.uplo = C->prop.uplo;
            herk_dnx_dny<T>::do_all(sub_A, sub_C, mode, alpha, T{1});
        }
    }
}



template<typename T, typename I>
void herk_dnx_dny_tri<T>::finalize()
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
