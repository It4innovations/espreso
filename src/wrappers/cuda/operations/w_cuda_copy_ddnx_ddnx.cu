
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cuda_copy_ddnx_ddnx.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
__global__
static void copy_popullatedfirst(T * src, size_t ld_src, T * dst, size_t ld_dst)
{
    size_t ip = blockIdx.x;
    T * sub_src = src + ld_src * ip;
    T * sub_dst = dst + ld_dst * ip;
    for(size_t is = threadIdx.x; is <= ip; is += blockDim.x) {
        sub_dst[is] = sub_src[is];
    }
}



template<typename T>
__global__
static void copy_emptyfirst(T * src, size_t ld_src, T * dst, size_t ld_dst, size_t size_secdary)
{
    size_t ip = blockIdx.x;
    T * sub_src = src + ld_src * ip;
    T * sub_dst = dst + ld_dst * ip;
    size_t start = (ip / warpSize) * warpSize;
    for(size_t is = start + threadIdx.x; is < size_secdary; is += blockDim.x) {
        if(is >= sp) {
            sub_dst[is] = sub_src[is];
        }
    }
}



template<typename T>
void w_cuda_copy_ddnx_ddnx<T>::internal_setup()
{
    wss_tmp_perform = 0;
}



template<typename T>
void w_cuda_copy_ddnx_ddnx<T>::internal_perform(void * ws_tmp)
{
    if(uplo == 'L' || uplo == 'U') {
        if((uplo == 'L') == (M_src->order == 'R')) {
            copy_popullatedfirst<T><<<M_src->get_size_primary(),256>>>(M_src->vals, M_src->ld, M_dst->vals, M_dst->ld);
            CHECK(cudaPeekAtLastError());
        }
        else {
            copy_emptyfirst<T><<<M_src->get_size_primary(),256>>>(M_src->vals, M_src->ld, M_dst->vals, M_dst->ld, M_src->get_size_secdary());
            CHECK(cudaPeekAtLastError());
        }
    }
    else {
        CHECK(cudaMemcpy2DAsync(M_dst->vals, M_dst->ld * sizeof(T), M_src->vals, M_src->ld * sizeof(T), M_src->get_size_secdary() * sizeof(T), M_src->get_size_primary(), cudaMemcpyDefault, q->stream));
    }
}



#define INSTANTIATE_T(T) \
template class w_cuda_copy_ddnx_ddnx<T>;

    #define INSTANTIATE \
    /* INSTANTIATE_T(float) */ \
    INSTANTIATE_T(double) \
    /* INSTANTIATE_T(std::complex<float>) */ \
    /* INSTANTIATE_T(std::complex<double>) */

        INSTANTIATE

    #undef INSTANTIATE
#undef INSTANTIATE_T



}
}
}

#endif
