
#include "gpu/operations/trsm_dcsx_ddny_ddny.h"

#include "gpu/gpu_spblas.h"

#include "wrappers/cuda/operations/w_cusparse_trsm_dcsx_ddny_ddny.h"
// #include "wrappers/rocm/operations/w_rocsparse_trsm_dcsx_ddny_ddny.h"
// #include "wrappers/oneapi/operations/w_oneapisparse_trsm_dcsx_ddny_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
std::unique_ptr<trsm_dcsx_ddny_ddny<T,I>> trsm_dcsx_ddny_ddny<T,I>::make()
{
    // feel free to make this runtime ifs based on ecf or env
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<w_cusparse_trsm_dcsx_ddny_ddny<T,I>>();
    #endif
    // #ifdef ESPRESO_USE_WRAPPER_GPU_ROCM
    //     return std::make_unique<w_rocsparse_trsm_dcsx_ddny_ddny<T,I>>();
    // #endif
    // #ifdef ESPRESO_USE_WRAPPER_GPU_ONEAPI
    //     return std::make_unique<w_oneapisparse_trsm_dcsx_ddny_ddny<T,I>>();
    // #endif
    eslog::error("wrapper for trsm_dcsx_ddny_ddny not available\n");
}



template<typename T, typename I>
char trsm_dcsx_ddny_ddny<T,I>::get_native_place()
{
    return this->internal_get_native_place();
}



template<typename T, typename I>
void trsm_dcsx_ddny_ddny<T,I>::set_handles(gpu::mgm::queue q_, gpu::spblas::handle spblas_handle_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    spblas_handle = spblas_handle_;

    called_set_handles = true;
}



template<typename T, typename I>
void trsm_dcsx_ddny_ddny<T,I>::set_matrix_A(MatrixCsxView_new<T,I> A_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(called_set_A) eslog::error("forbidden to re-set matrix A\n");

    A = A_;

    internal_set_matrix_A();

    called_set_A = true;
}



template<typename T, typename I>
void trsm_dcsx_ddny_ddny<T,I>::set_matrix_X(MatrixDenseView_new<T> X_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(called_set_X && !MatrixDenseView_new<T>::are_interchangable(X, X_)) eslog::error("invalid replacement for matrix X\n");

    X = X_;

    internal_set_matrix_X();

    called_set_X = true;
}



template<typename T, typename I>
void trsm_dcsx_ddny_ddny<T,I>::set_matrix_B(MatrixDenseView_new<T> B_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(called_set_B && !MatrixDenseView_new<T>::are_interchangable(B, B_)) eslog::error("invalid replacement for matrix B\n");

    B = B_;

    internal_set_matrix_B();

    called_set_B = true;
}



template<typename T, typename I>
void trsm_dcsx_ddny_ddny<T,I>::setup()
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(!called_set_A) eslog::error("matrix A is not set\n");
    if(!called_set_B) eslog::error("matrix B is not set\n");
    if(!called_set_C) eslog::error("matrix C is not set\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(A.nrows != A.ncols) eslog::error("matrix A is not square\n");
    if(X.nrows != B.nrows || X.ncols != B.ncols) eslog::error("X and B matrix sizes dont match\n");
    if(X.order != B.order) eslog::error("X and B order does not match\n");
    if(A.nrows != B.nrows) eslog::error("incompatible matrices\n");

    place = ((X.vals == B.vals) ? 'I' : 'O');

    this->internal_setup();

    called_setup = true;
}



template<typename T, typename I>
size_t trsm_dcsx_ddny_ddny<T,I>::get_wss_internal()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_internal;
}



template<typename T, typename I>
size_t trsm_dcsx_ddny_ddny<T,I>::get_wss_persistent()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_persistent;
}



template<typename T, typename I>
size_t trsm_dcsx_ddny_ddny<T,I>::get_wss_tmp_preprocess()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_preprocess;
}



template<typename T, typename I>
size_t trsm_dcsx_ddny_ddny<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void trsm_dcsx_ddny_ddny<T,I>::set_ws_persistent(void * ws_persistent_)
{
    if(ws_persistent_ == nullptr && wss_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_persistent = ws_persistent_;
}



template<typename T, typename I>
void trsm_dcsx_ddny_ddny<T,I>::preprocess_submit(void * ws_tmp)
{
    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    this->internal_preprocess(ws_tmp);

    called_preprocess = true;
}



template<typename T, typename I>
void trsm_dcsx_ddny_ddny<T,I>::update_submit()
{
    if(!called_preprocess) eslog::error("preprocess has not been called\n");

    this->internal_update();
}



template<typename T, typename I>
void trsm_dcsx_ddny_ddny<T,I>::perform_submit(void * ws_tmp)
{
    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    this->internal_perform(ws_tmp);
}



template<typename T, typename I>
void trsm_dcsx_ddny_ddny<T,I>::submit_all(gpu::mgm::queue q, gpu::spblas::handle spblas_handle, MatrixCsxView_new<T,I> A, MatrixDenseView_new<T> x, MatrixDenseView_new<T> B, Allocator_new * ator_gpu)
{
    trsm_dcsx_ddny_ddny<T,I> instance;
    instance.set_handles(q, spblas_handle);
    instance.set_matrix_A(A);
    instance.set_matrix_X(X);
    instance.set_matrix_B(B);
    instance.setup();
    size_t wss_persistent = instance.get_wss_persistent();
    size_t wss_tmp_preprocess = instance.get_wss_tmp_preprocess();
    size_t wss_tmp_perform = instance.get_wss_tmp_perform();
    size_t wss_tmp = std::max(wss_tmp_preprocess, wss_tmp_perform);
    void * ws_persistent = ator_gpu->alloc(wss_persistent);
    void * ws_tmp = ator_gpu->alloc(wss_tmp);
    instance.set_ws_persistent(ws_persistent);
    instance.preprocess_submit(ws_tmp);
    instance.perform_submit(ws_tmp);
    gpu::mgm::submit_host_function(q, [ator_gpu,ws_persistent,ws_tmp](){
        ator_gpu->free(ws_persistent);
        ator_gpu->free(ws_tmp);
    });
}



#define INSTANTIATE_T_I(T,I) \
template class trsm_dcsx_ddny_ddny<T,I>;

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
