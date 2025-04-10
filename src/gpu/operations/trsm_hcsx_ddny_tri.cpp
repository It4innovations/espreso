
#include "gpu/operations/trsm_hcsx_ddny_tri.h"

#include "basis/utilities/stacktimer.h"
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
void trsm_hcsx_ddny_tri<T,I>::set_matrix_h_L(MatrixCsxView_new<T,I> * h_L_)
{
    if(h_L != nullptr) eslog::error("matrix h_L is already set\n");
    if(h_L_ == nullptr) eslog::error("h_L cannot be nullptr\n");

    h_L = h_L_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::set_matrix_d_X(MatrixDenseView_new<T> * d_X_)
{
    if(d_X != nullptr) eslog::error("matrix d_X is already set\n");
    if(d_X_ == nullptr) eslog::error("d_X cannot be nullptr\n");

    d_X = d_X_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::set_X_pattern(MatrixCsxView_new<T,I> * h_X_pattern_)
{
    if(h_X_pattern != nullptr) eslog::error("X pattern is already set\n");
    if(h_X_pattern_ == nullptr) eslog::error("X pattern cannot be nullptr\n");

    h_X_pattern = h_X_pattern_;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::setup()
{
    stacktimer::push("trsm_hcsx_ddny_tri::setup");

    if(!called_set_handles) eslog::error("handles are not set\n");
    if(h_L == nullptr) eslog::error("matrix L is not set\n");
    if(d_X == nullptr) eslog::error("matrix X is not set\n");
    if(h_X_pattern == nullptr) eslog::error("X pattern is not set\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(h_L->nrows != h_L->ncols) eslog::error("matrix L is not square\n");
    if(h_L->nrows != d_X->nrows) eslog::error("incompatible matrices\n");
    if(h_L->prop.uplo != 'L') eslog::error("L has to have uplo=L\n");

    ator_ws_persistent = std::make_unique<AllocatorArena_new>(false, true, gpu::mgm::get_natural_pitch_align());
    ator_ws_tmp_linear = std::make_unique<AllocatorArena_new>(false, true, gpu::mgm::get_natural_pitch_align());
    ator_ws_tmp_overlap = std::make_unique<AllocatorSinglePointer_new>(false, true, gpu::mgm::get_natural_pitch_align());

    if(cfg.strategy == 'R') {
        typename trsm_hcsx_ddny_tri_splitrhs<T,I>::config op_config;
        op_config.partition.algorithm = cfg.partition.algorithm;
        op_config.partition.parameter = cfg.partition.parameter;
        op_config.factor_order_sp = cfg.splitrhs.factor_order_sp;
        op_config.factor_order_dn = cfg.splitrhs.factor_order_dn;
        op_config.spdn_criteria = cfg.splitrhs.spdn_criteria;
        op_config.spdn_param = cfg.splitrhs.spdn_param;
        trsm_splitrhs = std::make_unique<trsm_hcsx_ddny_tri_splitrhs<T,I>>();
        trsm_splitrhs->set_config(op_config);
        trsm_splitrhs->set_handles(q, handle_spblas, handle_dnblas);
        trsm_splitrhs->set_matrix_h_L(h_L);
        trsm_splitrhs->set_matrix_d_X(d_X);
        trsm_splitrhs->set_h_X_pattern(h_X_pattern);
        trsm_splitrhs->setup();
        wss_internal += trsm_splitrhs->get_wss_internal();
        wss_persistent += utils::round_up(trsm_splitrhs->get_wss_persistent(), ator_ws_persistent->get_align());
        wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, trsm_splitrhs->get_wss_tmp_preprocess());
        wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, trsm_splitrhs->get_wss_tmp_perform());
    }
    if(cfg.strategy == 'F') {
        typename trsm_hcsx_ddny_tri_splitfactor<T,I>::config op_config;
        op_config.partition.algorithm = cfg.partition.algorithm;
        op_config.partition.parameter = cfg.partition.parameter;
        op_config.trsm_factor_spdn = cfg.splitfactor.trsm_factor_spdn;
        op_config.trsm_factor_order = cfg.splitfactor.trsm_factor_order;
        op_config.gemm_factor_prune = cfg.splitfactor.gemm_factor_prune;
        op_config.gemm_factor_order_sp = cfg.splitfactor.gemm_factor_order_sp;
        op_config.gemm_factor_order_dn = cfg.splitfactor.gemm_factor_order_dn;
        op_config.gemm_spdn_criteria = cfg.splitfactor.gemm_spdn_criteria;
        op_config.gemm_spdn_param = cfg.splitfactor.gemm_spdn_param;
        trsm_splitfactor = std::make_unique<trsm_hcsx_ddny_tri_splitfactor<T,I>>();
        trsm_splitfactor->set_config(op_config);
        trsm_splitfactor->set_handles(q, handle_spblas, handle_dnblas);
        trsm_splitfactor->set_matrix_h_L(h_L);
        trsm_splitfactor->set_matrix_d_X(d_X);
        trsm_splitfactor->set_h_X_pattern(h_X_pattern);
        trsm_splitfactor->setup();
        wss_internal += trsm_splitfactor->get_wss_internal();
        wss_persistent += utils::round_up(trsm_splitfactor->get_wss_persistent(), ator_ws_persistent->get_align());
        wss_tmp_preprocess_overlap = std::max(wss_tmp_preprocess_overlap, trsm_splitfactor->get_wss_tmp_preprocess());
        wss_tmp_perform_overlap = std::max(wss_tmp_perform_overlap, trsm_splitfactor->get_wss_tmp_perform());
    }

    wss_tmp_preprocess_linear = utils::round_up(wss_tmp_preprocess_linear, ator_ws_tmp_linear->get_align());
    wss_tmp_perform_linear = utils::round_up(wss_tmp_perform_linear, ator_ws_tmp_linear->get_align());

    wss_tmp_preprocess = wss_tmp_preprocess_linear + wss_tmp_preprocess_overlap;
    wss_tmp_perform = wss_tmp_perform_linear + wss_tmp_perform_overlap;

    stacktimer::info("wss_internal       %zu", wss_internal);
    stacktimer::info("wss_persistent     %zu", wss_persistent);
    stacktimer::info("wss_tmp_preprocess %zu", wss_tmp_preprocess);
    stacktimer::info("wss_tmp_perform    %zu", wss_tmp_perform);

    stacktimer::pop();

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
    stacktimer::push("trsm_hcsx_ddny_tri::preprocess_submit");

    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    ator_ws_persistent->set(ws_persistent, wss_persistent);

    ator_ws_tmp_linear->set(ws_tmp, wss_tmp_preprocess_linear);
    ator_ws_tmp_overlap->set((char*)ws_tmp + wss_tmp_preprocess_linear, wss_tmp_preprocess_overlap);

    if(cfg.strategy == 'R') {
        trsm_splitrhs->set_ws_persistent(ator_ws_persistent->alloc(trsm_splitrhs->get_wss_persistent()));
        trsm_splitrhs->preprocess_submit(ator_ws_tmp_overlap->alloc(trsm_splitrhs->get_wss_tmp_preprocess()));
    }
    if(cfg.strategy == 'F') {
        trsm_splitfactor->set_ws_persistent(ator_ws_persistent->alloc(trsm_splitfactor->get_wss_persistent()));
        trsm_splitfactor->preprocess_submit(ator_ws_tmp_overlap->alloc(trsm_splitfactor->get_wss_tmp_preprocess()));
    }
    
    ator_ws_tmp_linear->unset();
    ator_ws_tmp_overlap->unset();

    stacktimer::pop();

    called_preprocess = true;
}



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::perform_submit(void * ws_tmp)
{
    stacktimer::push("trsm_hcsx_ddny_tri::perform_submit");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    ator_ws_tmp_linear->set(ws_tmp, wss_tmp_perform_linear);
    ator_ws_tmp_overlap->set((char*)ws_tmp + wss_tmp_perform_linear, wss_tmp_perform_overlap);

    if(cfg.strategy == 'R') {
        trsm_splitrhs->perform_submit(ator_ws_tmp_overlap->alloc(trsm_splitrhs->get_wss_tmp_perform()));
    }
    if(cfg.strategy == 'F') {
        trsm_splitfactor->perform_submit(ator_ws_tmp_overlap->alloc(trsm_splitfactor->get_wss_tmp_perform()));
    }

    ator_ws_tmp_linear->unset();
    ator_ws_tmp_overlap->unset();

    stacktimer::pop();
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
