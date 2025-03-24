
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cuda_convert_ddnx_ddny.h"

#include "wrappers/cuda/common_cuda_mgm.h"
#include "wrappers/cuda/common_cublas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
void w_cuda_convert_ddnx_ddny<T>::internal_setup()
{
    wss_tmp_perform = 0;
}



template<typename T>
void w_cuda_convert_ddnx_ddny<T>::internal_perform(void * ws_tmp)
{
    using U = cpp_to_cuda_type_t<T>;
    if(M_src->order == M_dst->order) {
        CHECK(cudaMemcpy2DAsync(M_dst->vals, M_dst->ld * sizeof(T), M_src->vals, M_src->ld * sizeof(T), M_src->get_size_secdary() * sizeof(T), M_src->get_size_primary(), cudaMemcpyDefault, q->stream));
    }
    else {
        T one = T{1};
        T zero = T{0};
        if constexpr(std::is_same_v<T,float>)                CHECK(cublasSgeam(handle_dnblas->h, CUBLAS_OP_T, CUBLAS_OP_N, M_src->get_size_secdary(), M_src->get_size_primary(), (U*)&one, (U*)M_src->vals, M_src->ld, (U*)&zero, nullptr, M_src->ld, (U*)M_dst->vals, M_dst->ld));
        if constexpr(std::is_same_v<T,double>)               CHECK(cublasDgeam(handle_dnblas->h, CUBLAS_OP_T, CUBLAS_OP_N, M_src->get_size_secdary(), M_src->get_size_primary(), (U*)&one, (U*)M_src->vals, M_src->ld, (U*)&zero, nullptr, M_src->ld, (U*)M_dst->vals, M_dst->ld));
        if constexpr(std::is_same_v<T,std::complex<float>>)  CHECK(cublasCgeam(handle_dnblas->h, CUBLAS_OP_T, CUBLAS_OP_N, M_src->get_size_secdary(), M_src->get_size_primary(), (U*)&one, (U*)M_src->vals, M_src->ld, (U*)&zero, nullptr, M_src->ld, (U*)M_dst->vals, M_dst->ld));
        if constexpr(std::is_same_v<T,std::complex<double>>) CHECK(cublasZgeam(handle_dnblas->h, CUBLAS_OP_T, CUBLAS_OP_N, M_src->get_size_secdary(), M_src->get_size_primary(), (U*)&one, (U*)M_src->vals, M_src->ld, (U*)&zero, nullptr, M_src->ld, (U*)M_dst->vals, M_dst->ld));
    }
}



#define INSTANTIATE_T(T) \
template class w_cuda_convert_ddnx_ddny<T>;

    #define INSTANTIATE \
    /* INSTANTIATE_T(float) */ \
    INSTANTIATE_T(double) \
    /* INSTANTIATE_T(std::complex<float>) */ \
    INSTANTIATE_T(std::complex<double>)

        INSTANTIATE

    #undef INSTANTIATE
#undef INSTANTIATE_T



}
}
}

#endif
