
#ifdef HAVE_CUDA

#include "wrappers/cuda/operations/w_cusparse_convert_dcsx_ddny.h"

#include "wrappers/cuda/common_cusparse.h"



namespace espreso {
namespace gpu {
namespace operations {



struct w_cusparse_convert_dcsx_ddny_data
{
    cusparseSpMatDescr_t descr_M_src;
    cusparseDnMatDescr_t descr_M_dst;
    cusparseSparseToDenseAlg_t alg;
};



template<typename T, typename I>
w_cusparse_convert_dcsx_ddny<T,I>::w_cusparse_convert_dcsx_ddny()
{
    data = std::make_unique<w_cusparse_convert_dcsx_ddny_data>();
}



template<typename T, typename I>
w_cusparse_convert_dcsx_ddny<T,I>::~w_cusparse_convert_dcsx_ddny()
{
    data.reset();
}



template<typename T, typename I>
void w_cusparse_convert_dcsx_ddny<T,I>::internal_setup()
{
    T * dummyptr_T = (T*)(sizeof(T));
    I * dummyptr_I = (I*)(sizeof(I));
    if(M_src->order == 'R') CHECK(cusparseCreateCsr(&data->descr_M_src, M_src->nrows, M_src->ncols, M_src->nnz, dummyptr_I, dummyptr_I, dummyptr_T, cusparse_index_type<I>(), cusparse_index_type<I>(), CUSPARSE_INDEX_BASE_ZERO, cusparse_data_type<T>()));
    if(M_src->order == 'C') CHECK(cusparseCreateCsc(&data->descr_M_src, M_src->nrows, M_src->ncols, M_src->nnz, dummyptr_I, dummyptr_I, dummyptr_T, cusparse_index_type<I>(), cusparse_index_type<I>(), CUSPARSE_INDEX_BASE_ZERO, cusparse_data_type<T>()));

    CHECK(cusparseCreateDnMat(&data->descr_M_dst, M_dst->nrows, M_dst->ncols, M_dst->ld, dummyptr_T, cusparse_data_type<T>(), cusparse_order(M_dst->order)));

    data->alg = CUSPARSE_SPARSETODENSE_ALG_DEFAULT;

    wss_internal = 0; // TODO: check
    wss_persistent = 0;
    wss_tmp_preprocess = 0;
    CHECK(cusparseSparseToDense_bufferSize(handle_spblas->h, data->descr_M_src, data->descr_M_dst, data->alg, &wss_tmp_perform));
}



template<typename T, typename I>
void w_cusparse_convert_dcsx_ddny<T,I>::internal_preprocess(void * /*ws_tmp*/)
{
    // no-op
}



template<typename T, typename I>
void w_cusparse_convert_dcsx_ddny<T,I>::internal_perform(void * ws_tmp)
{
    if(M_src->order == 'R') CHECK(cusparseCsrSetPointers(data->descr_M_src, M_src->ptrs, M_src->idxs, M_src->vals));
    if(M_src->order == 'C') CHECK(cusparseCscSetPointers(data->descr_M_src, M_src->ptrs, M_src->idxs, M_src->vals));
    CHECK(cusparseDnMatSetValues(data->descr_M_dst, M_dst->vals));

    CHECK(cusparseSparseToDense(handle_spblas->h, data->descr_M_src, data->descr_M_dst, data->alg, ws_tmp));
}



#define INSTANTIATE_T_I(T,I) \
template class w_cusparse_convert_dcsx_ddny<T,I>;

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

#endif
