
#include "gpu/operations/trsm_hcsx_ddny_tri.h"

#include "config/ecf/operations/gpu_trsm_hcsx_ddny_tria.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "basis/utilities/stacktimer.h"
#include "math/operations/pivots_trails_csx.h"



namespace espreso {
namespace gpu {
namespace operations {



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
    if(!h_L->ator->is_data_accessible_cpu()) eslog::error("matrix h_L must be cpu-accessible\n");
    if(!d_X->ator->is_data_accessible_gpu()) eslog::error("matrix d_X must be gpu-accessible\n");
    if(!h_X_pattern->ator->is_data_accessible_cpu()) eslog::error("matrix h_X_pattern must be cpu-accessible\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(h_L->nrows != h_L->ncols) eslog::error("matrix L is not square\n");
    if(h_L->nrows != d_X->nrows) eslog::error("incompatible matrices\n");
    if(h_L->prop.uplo != 'L') eslog::error("L has to have uplo=L\n");

    setup_config();

    ator_ws_persistent = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());
    ator_ws_tmp_linear = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());
    ator_ws_tmp_overlap = std::make_unique<AllocatorSinglePointer_new>(AllocatorGPU_new::get_singleton());

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



template<typename T, typename I>
void trsm_hcsx_ddny_tri<T,I>::setup_config()
{
    using ecf_config = GpuTrsmHcsxDdnyTriaConfig;
    const ecf_config & ecf = info::ecf->operations.gpu_trsm_hcsx_ddny_tria;

    switch(ecf.strategy) {
        case ecf_config::TRSM_TRIA_STRATEGY::AUTO:         cfg.strategy = 'F'; break;
        case ecf_config::TRSM_TRIA_STRATEGY::SPLIT_RHS:    cfg.strategy = 'R'; break;
        case ecf_config::TRSM_TRIA_STRATEGY::SPLIT_FACTOR: cfg.strategy = 'F'; break;
    }

    {
        switch(ecf.partition.algorithm) {
            case ecf_config::PARTITION_ALGORITHM::AUTO:         cfg.partition.algorithm = 'U'; break;
            case ecf_config::PARTITION_ALGORITHM::UNIFORM:      cfg.partition.algorithm = 'U'; break;
            case ecf_config::PARTITION_ALGORITHM::MINIMUM_WORK: cfg.partition.algorithm = 'M'; break;
        }

        char partition_strategy = '_';
        switch(ecf.partition.strategy) {
            case ecf_config::PARTITION_STRATEGY::AUTO: {
                if(info::mesh->dimension == 2 && cfg.strategy == 'F') partition_strategy = 'S';
                if(info::mesh->dimension == 2 && cfg.strategy == 'R') partition_strategy = 'C';
                if(info::mesh->dimension == 3 && cfg.strategy == 'F') partition_strategy = 'S';
                if(info::mesh->dimension == 3 && cfg.strategy == 'R') partition_strategy = 'S';
                break;
            }
            case ecf_config::PARTITION_STRATEGY::CHUNK_SIZE:  partition_strategy = 'S'; break;
            case ecf_config::PARTITION_STRATEGY::CHUNK_COUNT: partition_strategy = 'C'; break;
        }

        int chunk_size = ecf.partition.chunk_size;
        if(chunk_size == 0) {
            if(info::mesh->dimension == 2 && cfg.strategy == 'F') chunk_size = 1000;
            if(info::mesh->dimension == 2 && cfg.strategy == 'R') chunk_size = 1000; // not tested
            if(info::mesh->dimension == 3 && cfg.strategy == 'F') chunk_size = 500;
            if(info::mesh->dimension == 3 && cfg.strategy == 'R') chunk_size = 1000;
        }

        int chunk_count = utils::replace_if_zero(ecf.partition.chunk_count, 1); // not tested

        if(partition_strategy == 'S') cfg.partition.parameter = -chunk_size;
        if(partition_strategy == 'C') cfg.partition.parameter = chunk_count;
    }

    {
        switch(ecf.split_rhs_config.factor_order_sp) {
            case ecf_config::MATRIX_ORDER::AUTO:      cfg.splitrhs.factor_order_sp = 'R'; break;
            case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.splitrhs.factor_order_sp = 'R'; break;
            case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.splitrhs.factor_order_sp = 'C'; break;
        }

        switch(ecf.split_rhs_config.factor_order_dn) {
            case ecf_config::MATRIX_ORDER::AUTO:      cfg.splitrhs.factor_order_dn = 'C'; break;
            case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.splitrhs.factor_order_dn = 'R'; break;
            case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.splitrhs.factor_order_dn = 'C'; break;
        }

        switch(ecf.split_rhs_config.spdn_criteria) {
            case ecf_config::SPDN_CRITERIA::AUTO:                    cfg.splitrhs.spdn_criteria = 'S'; break;
            case ecf_config::SPDN_CRITERIA::SPARSE_ONLY:             cfg.splitrhs.spdn_criteria = 'S'; break;
            case ecf_config::SPDN_CRITERIA::DENSE_ONLY:              cfg.splitrhs.spdn_criteria = 'D'; break;
            case ecf_config::SPDN_CRITERIA::FRACTION_OF_NUM_CHUNKS:  cfg.splitrhs.spdn_criteria = 'C'; break;
            case ecf_config::SPDN_CRITERIA::FRACTION_OF_FACTOR_SIZE: cfg.splitrhs.spdn_criteria = 'Z'; break;
            case ecf_config::SPDN_CRITERIA::FACTOR_DENSITY:          cfg.splitrhs.spdn_criteria = 'T'; break;
        }
        
        double spdn_param_frac_of_num_chunks = utils::replace_if_zero(ecf.split_rhs_config.spdn_param_frac_of_num_chunks, 0.7); // not tested
        double spdn_param_frac_of_factor_size = utils::replace_if_zero(ecf.split_rhs_config.spdn_param_frac_of_factor_size, 0.7); // not tested
        double spdn_param_factor_density = utils::replace_if_zero(ecf.split_rhs_config.spdn_param_factor_density, 0.1); // not tested

        cfg.splitrhs.spdn_param = 0;
        if(cfg.splitrhs.spdn_criteria == 'C') cfg.splitrhs.spdn_param = spdn_param_frac_of_num_chunks;
        if(cfg.splitrhs.spdn_criteria == 'Z') cfg.splitrhs.spdn_param = spdn_param_frac_of_factor_size;
        if(cfg.splitrhs.spdn_criteria == 'T') cfg.splitrhs.spdn_param = spdn_param_factor_density;
    }

    {
        switch(ecf.split_factor_config.trsm_factor_spdn) {
            case ecf_config::SPDN::AUTO:
                cfg.splitfactor.trsm_factor_spdn = ((info::mesh->dimension == 2) ? 'S' : 'D');
                break;
            case ecf_config::SPDN::SPARSE: cfg.splitfactor.trsm_factor_spdn = 'S'; break;
            case ecf_config::SPDN::DENSE:  cfg.splitfactor.trsm_factor_spdn = 'D'; break;
        }

        switch(ecf.split_factor_config.trsm_factor_order) {
            case ecf_config::MATRIX_ORDER::AUTO:      cfg.splitfactor.trsm_factor_order = 'R'; break;
            case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.splitfactor.trsm_factor_order = 'R'; break;
            case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.splitfactor.trsm_factor_order = 'C'; break;
        }

        switch(ecf.split_factor_config.gemm_factor_pruning) {
            case ecf_config::PRUNING_STRATEGY::AUTO:          cfg.splitfactor.gemm_factor_prune = 'R'; break;
            case ecf_config::PRUNING_STRATEGY::NO_PRUNING:    cfg.splitfactor.gemm_factor_prune = 'N'; break;
            case ecf_config::PRUNING_STRATEGY::ROWS_ONLY:     cfg.splitfactor.gemm_factor_prune = 'R'; break;
            case ecf_config::PRUNING_STRATEGY::COLS_ONLY:     cfg.splitfactor.gemm_factor_prune = 'C'; break;
            case ecf_config::PRUNING_STRATEGY::ROWS_AND_COLS: cfg.splitfactor.gemm_factor_prune = 'A'; break;
        }

        switch(ecf.split_factor_config.gemm_factor_order_sp) {
            case ecf_config::MATRIX_ORDER::AUTO:      cfg.splitfactor.gemm_factor_order_sp = 'R'; break;
            case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.splitfactor.gemm_factor_order_sp = 'R'; break;
            case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.splitfactor.gemm_factor_order_sp = 'C'; break;
        }

        switch(ecf.split_factor_config.gemm_factor_order_dn) {
            case ecf_config::MATRIX_ORDER::AUTO:      cfg.splitfactor.gemm_factor_order_dn = 'R'; break;
            case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.splitfactor.gemm_factor_order_dn = 'R'; break;
            case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.splitfactor.gemm_factor_order_dn = 'C'; break;
        }

        switch(ecf.split_factor_config.gemm_spdn_criteria) {
            case ecf_config::SPDN_CRITERIA::AUTO:
                cfg.splitfactor.gemm_spdn_criteria = ((info::mesh->dimension == 3 && cfg.splitfactor.gemm_factor_prune == 'R') ? 'D' : 'S');
                break;
            case ecf_config::SPDN_CRITERIA::SPARSE_ONLY:             cfg.splitfactor.gemm_spdn_criteria = 'S'; break;
            case ecf_config::SPDN_CRITERIA::DENSE_ONLY:              cfg.splitfactor.gemm_spdn_criteria = 'D'; break;
            case ecf_config::SPDN_CRITERIA::FRACTION_OF_NUM_CHUNKS:  cfg.splitfactor.gemm_spdn_criteria = 'C'; break;
            case ecf_config::SPDN_CRITERIA::FRACTION_OF_FACTOR_SIZE: cfg.splitfactor.gemm_spdn_criteria = 'Z'; break;
            case ecf_config::SPDN_CRITERIA::FACTOR_DENSITY:          cfg.splitfactor.gemm_spdn_criteria = 'T'; break;
        }
        
        double spdn_param_frac_of_num_chunks = utils::replace_if_zero(ecf.split_factor_config.spdn_param_frac_of_num_chunks, 0.7); // not tested
        double spdn_param_frac_of_factor_size = utils::replace_if_zero(ecf.split_factor_config.spdn_param_frac_of_factor_size, 0.7); // not tested
        double spdn_param_factor_density = utils::replace_if_zero(ecf.split_factor_config.spdn_param_factor_density, 0.1); // not tested

        cfg.splitfactor.gemm_spdn_param = 0;
        if(cfg.splitfactor.gemm_spdn_criteria == 'C') cfg.splitfactor.gemm_spdn_param = spdn_param_frac_of_num_chunks;
        if(cfg.splitfactor.gemm_spdn_criteria == 'Z') cfg.splitfactor.gemm_spdn_param = spdn_param_frac_of_factor_size;
        if(cfg.splitfactor.gemm_spdn_criteria == 'T') cfg.splitfactor.gemm_spdn_param = spdn_param_factor_density;
    }
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
