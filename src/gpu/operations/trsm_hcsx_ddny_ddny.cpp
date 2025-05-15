
#include "gpu/operations/trsm_hcsx_ddny_ddny.h"

#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
void trsm_hcsx_ddny_ddny<T,I>::set_config(char spdn_A_)
{
    if(called_set_config) eslog::error("config is already set\n");

    spdn_A = spdn_A_;

    called_set_config = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_ddny<T,I>::set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_, gpu::dnblas::handle handle_dnblas_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_spblas = handle_spblas_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_ddny<T,I>::set_matrix_h_A(MatrixCsxView_new<T,I> * h_A_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(h_A != nullptr) eslog::error("matrix h_A is already set\n");
    if(h_A_ == nullptr) eslog::error("h_A cannot be nullptr\n");

    h_A = h_A_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_ddny<T,I>::set_matrix_d_X(MatrixDenseView_new<T> * d_X_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(d_X != nullptr) eslog::error("matrix d_X is already set\n");
    if(d_X_ == nullptr) eslog::error("d_X cannot be nullptr\n");

    d_X = d_X_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_ddny<T,I>::set_matrix_d_B(MatrixDenseView_new<T> * d_B_)
{
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(d_B != nullptr) eslog::error("matrix d_B is already set\n");
    if(d_B_ == nullptr) eslog::error("d_B cannot be nullptr\n");

    d_B = d_B_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_ddny<T,I>::setup()
{
    stacktimer::push("trsm_hcsx_ddny_ddny::setup");

    if(!called_set_config) eslog::error("config is not set\n");
    if(!called_set_handles) eslog::error("handles are not set\n");
    if(h_A == nullptr) eslog::error("matrix h_A is not set\n");
    if(d_X == nullptr) eslog::error("matrix d_B is not set\n");
    if(d_B == nullptr) eslog::error("matrix d_C is not set\n");
    if(!h_A->ator->is_data_accessible_cpu()) eslog::error("matrix h_A must be cpu-accessible\n");
    if(!d_X->ator->is_data_accessible_gpu()) eslog::error("matrix d_X must be gpu-accessible\n");
    if(!d_B->ator->is_data_accessible_gpu()) eslog::error("matrix d_B must be gpu-accessible\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(h_A->nrows != h_A->ncols) eslog::error("matrix A is not square\n");
    if(d_X->nrows != d_B->nrows || d_X->ncols != d_B->ncols) eslog::error("X and B matrix sizes dont match\n");
    if(d_X->order != d_B->order) eslog::error("X and B order does not match\n");
    if(h_A->nrows != d_B->nrows) eslog::error("incompatible matrices\n");
    if(h_A->prop.uplo != 'L' && h_A->prop.uplo != 'U') eslog::error("matrix h_A has wrong uplo\n");

    ator_ws_persistent = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());

    ator_ws_tmp_linear = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());
    ator_ws_tmp_overlap = std::make_unique<AllocatorSinglePointer_new>(AllocatorGPU_new::get_singleton());

    d_A_sp.set(h_A->nrows, h_A->ncols, h_A->nnz, h_A->order, ator_ws_tmp_linear.get());
    d_A_sp.prop.uplo = h_A->prop.uplo;
    d_A_sp.prop.diag = h_A->prop.diag;
    wss_tmp_preprocess_linear += d_A_sp.get_memory_impact();
    wss_tmp_perform_linear += d_A_sp.get_memory_impact();

    if(spdn_A == 'D') {
        d_A_dn.set(d_A_sp.nrows, d_A_sp.ncols, d_A_sp.order, ator_ws_tmp_linear.get());
        d_A_dn.prop.uplo = h_A->prop.uplo;
        d_A_dn.prop.diag = h_A->prop.diag;
        wss_tmp_perform_linear += d_A_dn.get_memory_impact();
    }

    if(spdn_A == 'D') {
        op_d_sp2dn_A = convert_dcsx_ddny<T,I>::make();
        op_d_sp2dn_A->set_handles(q, handle_spblas);
        op_d_sp2dn_A->set_matrix_src(&d_A_sp);
        op_d_sp2dn_A->set_matrix_dst(&d_A_dn);
        op_d_sp2dn_A->setup();
        wss_internal += op_d_sp2dn_A->get_wss_internal();
        wss_persistent += utils::round_up(op_d_sp2dn_A->get_wss_persistent(), ator_ws_persistent->get_align());
        wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, op_d_sp2dn_A->get_wss_tmp_preprocess());
        wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_d_sp2dn_A->get_wss_tmp_perform());
    }

    if(spdn_A == 'S') {
        op_d_inner_trsm_sp = trsm_dcsx_ddny_ddny<T,I>::make();
        op_d_inner_trsm_sp->set_handles(q, handle_spblas);
        op_d_inner_trsm_sp->set_matrix_A(&d_A_sp);
        op_d_inner_trsm_sp->set_matrix_X(d_X);
        op_d_inner_trsm_sp->set_matrix_B(d_B);
        op_d_inner_trsm_sp->setup();
        wss_internal += op_d_inner_trsm_sp->get_wss_internal();
        wss_persistent += utils::round_up(op_d_inner_trsm_sp->get_wss_persistent(), ator_ws_persistent->get_align());
        wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, op_d_inner_trsm_sp->get_wss_tmp_preprocess());
        wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_d_inner_trsm_sp->get_wss_tmp_perform());
    }
    if(spdn_A == 'D') {
        op_d_inner_trsm_dn = trsm_ddnx_ddny<T>::make();
        op_d_inner_trsm_dn->set_handles(q, handle_dnblas);
        op_d_inner_trsm_dn->set_matrix_A(&d_A_dn);
        op_d_inner_trsm_dn->set_matrix_X(d_X);
        op_d_inner_trsm_dn->setup();
        wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, op_d_inner_trsm_dn->get_wss_tmp_perform());
    }

    wss_tmp_preprocess_linear = ((wss_tmp_preprocess_linear - 1) / ator_ws_tmp_linear->get_align() + 1) * ator_ws_tmp_linear->get_align();
    wss_tmp_perform_linear = ((wss_tmp_perform_linear - 1) / ator_ws_tmp_linear->get_align() + 1) * ator_ws_tmp_linear->get_align();

    wss_tmp_preprocess = wss_tmp_preprocess_linear + wss_tmp_preprocess_overlap;
    wss_tmp_perform = wss_tmp_perform_linear + wss_tmp_perform_overlap;

    // stacktimer::info("wss_internal       %zu", wss_internal);
    // stacktimer::info("wss_persistent     %zu", wss_persistent);
    // stacktimer::info("wss_tmp_preprocess %zu", wss_tmp_preprocess);
    // stacktimer::info("wss_tmp_perform    %zu", wss_tmp_perform);

    stacktimer::pop();

    called_setup = true;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_ddny<T,I>::get_wss_internal()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_internal;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_ddny<T,I>::get_wss_persistent()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_persistent;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_ddny<T,I>::get_wss_tmp_preprocess()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_preprocess;
}



template<typename T, typename I>
size_t trsm_hcsx_ddny_ddny<T,I>::get_wss_tmp_perform()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return wss_tmp_perform;
}



template<typename T, typename I>
void trsm_hcsx_ddny_ddny<T,I>::set_ws_persistent(void * ws_persistent_)
{
    if(ws_persistent_ == nullptr && wss_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_persistent = ws_persistent_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_ddny<T,I>::preprocess_submit(void * ws_tmp)
{
    stacktimer::push("trsm_hcsx_ddny_ddny::preprocess_submit");

    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    ator_ws_persistent->set(ws_persistent, wss_persistent);

    ator_ws_tmp_linear->set(ws_tmp, wss_tmp_preprocess_linear);
    ator_ws_tmp_overlap->set((char*)ws_tmp + wss_tmp_preprocess_linear, wss_tmp_preprocess_overlap);

    d_A_sp.alloc();

    gpu::mgm::copy_submit(q, *h_A, d_A_sp, true, false);

    if(spdn_A == 'D') {
        op_d_sp2dn_A->set_ws_persistent(ator_ws_persistent->alloc(op_d_sp2dn_A->get_wss_persistent()));
        op_d_sp2dn_A->preprocess_submit(ator_ws_tmp_overlap->alloc(op_d_sp2dn_A->get_wss_tmp_preprocess()));
    }

    if(spdn_A == 'S') {
        op_d_inner_trsm_sp->set_ws_persistent(ator_ws_persistent->alloc(op_d_inner_trsm_sp->get_wss_persistent()));
        op_d_inner_trsm_sp->preprocess_submit(ator_ws_tmp_overlap->alloc(op_d_inner_trsm_sp->get_wss_tmp_preprocess()));
    }

    d_A_sp.free();

    ator_ws_tmp_linear->unset();
    ator_ws_tmp_overlap->unset();

    stacktimer::pop();

    called_preprocess = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_ddny<T,I>::perform_submit(void * ws_tmp)
{
    stacktimer::push("trsm_hcsx_ddny_ddny::perform_submit");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    ator_ws_tmp_linear->set(ws_tmp, wss_tmp_perform_linear);
    ator_ws_tmp_overlap->set((char*)ws_tmp + wss_tmp_perform_linear, wss_tmp_perform_overlap);

    d_A_sp.alloc();
    if(spdn_A == 'D') {
        d_A_dn.alloc();
    }

    gpu::mgm::copy_submit(q, *h_A, d_A_sp, true, true);

    if(spdn_A == 'D') {
        op_d_sp2dn_A->perform_submit(ator_ws_tmp_overlap->alloc(op_d_sp2dn_A->get_wss_tmp_perform()));
    }

    if(spdn_A == 'S') {
        op_d_inner_trsm_sp->perform_submit(ator_ws_tmp_overlap->alloc(op_d_inner_trsm_sp->get_wss_tmp_perform()));
    }
    if(spdn_A == 'D') {
        op_d_inner_trsm_dn->perform_submit(ator_ws_tmp_overlap->alloc(op_d_inner_trsm_dn->get_wss_tmp_perform()));
    }

    d_A_sp.free();
    d_A_dn.free();

    ator_ws_tmp_linear->unset();
    ator_ws_tmp_overlap->unset();

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class trsm_hcsx_ddny_ddny<T,I>;

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
