
#ifdef HAVE_ROCM

#include "wrappers/rocm/operations/w_rocsparse_convert_dcsx_ddny.h"

#include "wrappers/rocm/common_rocsparse.h"



namespace espreso {
namespace gpu {
namespace operations {



struct w_rocsparse_convert_dcsx_ddny_data
{
    rocsparse_spmat_descr descr_M_src;
    rocsparse_dnmat_descr descr_M_dst;
    rocsparse_sparse_to_dense_alg alg;
};



template<typename T, typename I>
w_rocsparse_convert_dcsx_ddny<T,I>::w_rocsparse_convert_dcsx_ddny() {}



template<typename T, typename I>
w_rocsparse_convert_dcsx_ddny<T,I>::~w_rocsparse_convert_dcsx_ddny()
{
    if(data) {
        CHECK(rocsparse_destroy_spmat_descr(data->descr_M_src));
        CHECK(rocsparse_destroy_dnmat_descr(data->descr_M_dst));
    }
}



template<typename T, typename I>
void w_rocsparse_convert_dcsx_ddny<T,I>::internal_setup()
{
    data = std::make_unique<w_rocsparse_convert_dcsx_ddny_data>();

    T * dummyptr_T = (T*)(sizeof(T));
    I * dummyptr_I = (I*)(sizeof(I));
    if(M_src->order == 'R') CHECK(rocsparse_create_csr_descr(&data->descr_M_src, M_src->nrows, M_src->ncols, M_src->nnz, dummyptr_I, dummyptr_I, dummyptr_T, get_rocsparse_index_type<I>(), get_rocsparse_index_type<I>(), rocsparse_index_base_zero, get_rocsparse_data_type<T>()));
    if(M_src->order == 'C') CHECK(rocsparse_create_csc_descr(&data->descr_M_src, M_src->nrows, M_src->ncols, M_src->nnz, dummyptr_I, dummyptr_I, dummyptr_T, get_rocsparse_index_type<I>(), get_rocsparse_index_type<I>(), rocsparse_index_base_zero, get_rocsparse_data_type<T>()));

    CHECK(rocsparse_create_dnmat_descr(&data->descr_M_dst, M_dst->nrows, M_dst->ncols, M_dst->ld, dummyptr_T, get_rocsparse_data_type<T>(), get_rocsparse_order(M_dst->order)));

    data->alg = rocsparse_sparse_to_dense_alg_default;

    wss_internal = 0;
    wss_persistent = 0;
    wss_tmp_preprocess = 0;
    CHECK(rocsparse_sparse_to_dense(handle_spblas->h, data->descr_M_src, data->descr_M_dst, data->alg, &wss_tmp_perform, nullptr));
}



template<typename T, typename I>
void w_rocsparse_convert_dcsx_ddny<T,I>::internal_preprocess(void * /*ws_tmp*/)
{
    // no-op
}



template<typename T, typename I>
void w_rocsparse_convert_dcsx_ddny<T,I>::internal_perform(void * ws_tmp)
{
    if(M_src->order == 'R') CHECK(rocsparse_csr_set_pointers(data->descr_M_src, M_src->ptrs, M_src->idxs, M_src->vals));
    if(M_src->order == 'C') CHECK(rocsparse_csc_set_pointers(data->descr_M_src, M_src->ptrs, M_src->idxs, M_src->vals));
    CHECK(rocsparse_dnmat_set_values(data->descr_M_dst, M_dst->vals));

    CHECK(rocsparse_sparse_to_dense(handle_spblas->h, data->descr_M_src, data->descr_M_dst, data->alg, &wss_tmp_perform, ws_tmp));
}



#define INSTANTIATE_T_I(T,I) \
template class w_rocsparse_convert_dcsx_ddny<T,I>;

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
