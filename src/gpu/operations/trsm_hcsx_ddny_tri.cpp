
#include "gpu/operations/trsm_hcsx_ddny_tri.h"

#include "math/operations/pivots_trails_csx.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::set_config(config cfg_)
{
    cfg = cfg_;

    called_set_config = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::set_handles(gpu::mgm::queue q_, gpu::spblas::handle spblas_handle_, gpu::dnblas::handle dnblas_handle_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_spblas = spblas_handle_;
    handle_dnblas = dnblas_handle_;

    called_set_handles = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::set_matrix_L(MatrixCsxView_new<T,I> L)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(called_set_L) eslog::error("forbidden to re-set matrix L\n");

    L = L_;

    called_set_L = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::set_matrix_X(MatrixDenseView_new<T> X)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(called_set_X && !MatrixDenseView_new<T>::are_interchangable(X, X_)) eslog::error("invalid replacement for matrix X\n");

    X = X_;

    called_set_X = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::calc_X_pattern(MatrixCsxView<T,I> X_pattern_host)
{
    if(called_calc_X_pattern) eslog::error("X patern was already calculated\n");

    X_colpivots.set(X_pattern_host.ncols, AllocatorCPU_new::get_singleton());
    X_colpivots.alloc();
    pivots_trails_csx<T,I>::do_all(&X_pattern_host, &X_colpivots, 'C', 'P', 'B');

    X_rowtrails.set(X_pattern_host.nrows, AllocatorCPU_new::get_singleton());
    X_rowtrails.alloc();
    pivots_trails_csx<T,I>::do_all(&X_pattern_host, &X_rowtrails, 'R', 'T', 'F');

    for(size_t i = 1; i < X_colpivots.size; i++) {
        if(X_colpivots.vals[i-1] > X_colpivots.vals[i]) {
            eslog::error("X does not have lower triangular structure\n");
        }
    }
    for(size_t i = 1; i < X_rowtrails.size; i++) {
        if(X_rowtrails.vals[i-1] > X_rowtrails.vals[i]) {
            eslog::error("X does not have lower triangular structure\n");
        }
    }

    called_calc_X_pattern = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::setup()
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(!called_set_L) eslog::error("matrix L is not set\n");
    if(!called_set_X) eslog::error("matrix X is not set\n");
    if(!called_calc_X_pattern) eslog::error("X pattern has not been calculated\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(L.nrows != L.ncols) eslog::error("matrix L is not square\n");
    if(L.nrows != X.nrows) eslog::error("incompatible matrices\n");

    tri_partition_trsm partitioner;
    char partition_direction = '_';
    if(cfg.strategy == 'F') partition_direction = 'V';
    if(cfg.strategy == 'R') partition_direction = 'H';
    partitioner.set_config(cfg.partition.algorithm, partition_direction, cfg.partition.parameter);
    partitioner.set_system(X->nrows, X->ncols);
    partitioner.set_output_partition(&partition);
    partitioner.setup();
    num_chunks = partitioner.get_num_chunks();
    partition.set(num_chunks + 1, AllocatorCPU_new::get_singleton());
    partition.alloc();
    partitioner.perform();

    // TODO
    // wss_internal = 
    // wss_persistent = 
    // wss_tmp_preprocess = 
    // wss_tmp_perform = 

    called_setup = true;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_tri<T,I>::get_wss_internal()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_internal;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_tri<T,I>::get_wss_persistent()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_persistent;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_tri<T,I>::get_wss_tmp_preprocess()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_preprocess;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_tri<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::set_ws_persistent(void * ws_persistent_)
{
    if(ws_persistent_ == nullptr && wss_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_persistent = ws_persistent_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::preprocess_submit(void * ws_tmp)
{
    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    // TODO

    called_preprocess = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::update_submit()
{
    if(!called_preprocess) eslog::error("preprocess has not been called\n");

    // TODO
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::perform_submit(void * ws_tmp)
{
    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    // TODO

    // for all chunks perform chunk
}



#define INSTANTIATE_T_I(T,I) \
template class trsm_hcsx_ddny_tri<T,I>;

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
