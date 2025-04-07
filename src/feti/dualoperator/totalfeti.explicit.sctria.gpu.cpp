
#include "totalfeti.explicit.sctria.gpu.h"
#include "math/wrappers/math.blas.h"
#include "feti/common/applyB.h"
#include "basis/utilities/minmaxavg.h"
#include "my_timer.h"
#include "gpu/gpu_kernels.h"
#include "basis/utilities/stacktimer.h"
#include "esinfo/meshinfo.h"

#include <algorithm>



namespace espreso {

template<typename T>
static void set_by_env(T & var, const char * env_var)
{
    const char * str = std::getenv(env_var);
    if(str != nullptr) {
        if constexpr(std::is_same_v<T,char>) var = *str;
        if constexpr(std::is_same_v<T,int>) var = atoi(str);
        if constexpr(std::is_same_v<T,double>) var = atof(str);
        if constexpr(std::is_same_v<T,bool>) {
            var = (strcmp(str,"1") == 0
                || strcmp(str,"Y") == 0
                || strcmp(str,"Yes") == 0
                || strcmp(str,"yes") == 0
                || strcmp(str,"T") == 0
                || strcmp(str,"True") == 0
                || strcmp(str,"true") == 0
                || strcmp(str,"On") == 0
                || strcmp(str,"on") == 0);
        }
    }
}

static void replace_if_default(char & param, char deflt)
{
    if(param == '_') {
        param = deflt;
    }
}

static void replace_if_zero(int & param, int deflt)
{
    if(param == 0) {
        param = deflt;
    }
}

template<typename T, typename I>
static void replace_unset_configs(typename gpu::operations::sc_symm_hcsx_ddny_tria<T,I>::config & cfg_sc, typename TotalFETIExplicitScTriaGpu<T,I>::config & cfg_dualop)
{
    if(info::mesh->dimension == 2 && cfg_sc.cfg_herk.strategy == 'Q') replace_if_zero(cfg_sc.cfg_herk.partition_parameter, -500);
    if(info::mesh->dimension == 2 && cfg_sc.cfg_herk.strategy == 'T') replace_if_zero(cfg_sc.cfg_herk.partition_parameter, -200);
    if(info::mesh->dimension == 3 && cfg_sc.cfg_herk.strategy == 'Q') replace_if_zero(cfg_sc.cfg_herk.partition_parameter, -1000);
    if(info::mesh->dimension == 3 && cfg_sc.cfg_herk.strategy == 'T') replace_if_zero(cfg_sc.cfg_herk.partition_parameter, -1000);
    replace_if_default(cfg_sc.cfg_herk.partition_algorithm, 'U');

    if(info::mesh->dimension == 2 && cfg_sc.cfg_trsm.strategy == 'F') replace_if_zero(cfg_sc.cfg_trsm.partition.parameter, -1000);
    if(info::mesh->dimension == 2 && cfg_sc.cfg_trsm.strategy == 'R') replace_if_zero(cfg_sc.cfg_trsm.partition.parameter, 1);
    if(info::mesh->dimension == 3 && cfg_sc.cfg_trsm.strategy == 'F') replace_if_zero(cfg_sc.cfg_trsm.partition.parameter, -500);
    if(info::mesh->dimension == 3 && cfg_sc.cfg_trsm.strategy == 'R') replace_if_zero(cfg_sc.cfg_trsm.partition.parameter, -1000);
    replace_if_default(cfg_sc.cfg_trsm.partition.algorithm, 'U');

    replace_if_default(cfg_sc.cfg_trsm.splitrhs.spdn_criteria, 'S');
    replace_if_default(cfg_sc.cfg_trsm.splitfactor.trsm_factor_spdn, (info::mesh->dimension == 2) ? 'S' : 'D');
    replace_if_default(cfg_sc.cfg_trsm.splitfactor.gemm_factor_prune, 'R');
    replace_if_default(cfg_sc.cfg_trsm.splitfactor.gemm_spdn_criteria, (info::mesh->dimension == 2) ? 'S' : 'D');

    replace_if_default(cfg_sc.order_X, 'R');
    replace_if_default(cfg_sc.cfg_trsm.splitrhs.factor_order_sp, 'R');
    replace_if_default(cfg_sc.cfg_trsm.splitrhs.factor_order_dn, 'R');
    replace_if_default(cfg_sc.cfg_trsm.splitfactor.trsm_factor_order, 'R');
    replace_if_default(cfg_sc.cfg_trsm.splitfactor.gemm_factor_order_sp, 'R');
    replace_if_default(cfg_sc.cfg_trsm.splitfactor.gemm_factor_order_dn, 'R');

    replace_if_default(cfg_sc.order_L, 'C');

    replace_if_default(cfg_dualop.order_F, 'R');
    replace_if_default(cfg_dualop.mainloop_update_split, 'C');



    // replace_if_default(cfg_sc.cfg_trsm.strategy, '_');
    // replace_if_default(cfg_sc.cfg_herk.strategy, '_');
}

template<typename T, typename I>
static void setup_configs(typename gpu::operations::sc_symm_hcsx_ddny_tria<T,I>::config & cfg_sc, typename TotalFETIExplicitScTriaGpu<T,I>::config & cfg_dualop)
{
    cfg_sc.cfg_trsm.partition.parameter = 0;
    cfg_sc.cfg_herk.partition_parameter = 0;

    set_by_env(cfg_sc.order_X,                                   "ESPRESO_DUALOPSCTRIA_CONFIG_order_X");
    set_by_env(cfg_sc.order_L,                                   "ESPRESO_DUALOPSCTRIA_CONFIG_order_L");
    set_by_env(cfg_sc.cfg_trsm.strategy,                         "ESPRESO_DUALOPSCTRIA_CONFIG_trsm_strategy");
    set_by_env(cfg_sc.cfg_trsm.partition.algorithm,              "ESPRESO_DUALOPSCTRIA_CONFIG_trsm_partition_algorithm");
    set_by_env(cfg_sc.cfg_trsm.partition.parameter,              "ESPRESO_DUALOPSCTRIA_CONFIG_trsm_partition_parameter");
    set_by_env(cfg_sc.cfg_trsm.splitrhs.factor_order_sp,         "ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitrhs_factor_order_sp");
    set_by_env(cfg_sc.cfg_trsm.splitrhs.factor_order_dn,         "ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitrhs_factor_order_dn");
    set_by_env(cfg_sc.cfg_trsm.splitrhs.spdn_criteria,           "ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitrhs_spdn_criteria");
    set_by_env(cfg_sc.cfg_trsm.splitrhs.spdn_param,              "ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitrhs_spdn_param");
    set_by_env(cfg_sc.cfg_trsm.splitfactor.trsm_factor_spdn,     "ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitfactor_trsm_factor_spdn");
    set_by_env(cfg_sc.cfg_trsm.splitfactor.trsm_factor_order,    "ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitfactor_trsm_factor_order");
    set_by_env(cfg_sc.cfg_trsm.splitfactor.gemm_factor_order_sp, "ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitfactor_gemm_factor_order_sp");
    set_by_env(cfg_sc.cfg_trsm.splitfactor.gemm_factor_order_dn, "ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitfactor_gemm_factor_order_dn");
    set_by_env(cfg_sc.cfg_trsm.splitfactor.gemm_factor_prune,    "ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitfactor_gemm_factor_prune");
    set_by_env(cfg_sc.cfg_trsm.splitfactor.gemm_spdn_criteria,   "ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitfactor_gemm_spdn_criteria");
    set_by_env(cfg_sc.cfg_trsm.splitfactor.gemm_spdn_param,      "ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitfactor_gemm_spdn_param");
    set_by_env(cfg_sc.cfg_herk.strategy,                         "ESPRESO_DUALOPSCTRIA_CONFIG_herk_strategy");
    set_by_env(cfg_sc.cfg_herk.partition_algorithm,              "ESPRESO_DUALOPSCTRIA_CONFIG_herk_partition_algorithm");
    set_by_env(cfg_sc.cfg_herk.partition_parameter,              "ESPRESO_DUALOPSCTRIA_CONFIG_herk_partition_parameter");

    set_by_env(cfg_dualop.order_F,                        "ESPRESO_DUALOPSCTRIA_CONFIG_order_F");
    set_by_env(cfg_dualop.parallel_set,                   "ESPRESO_DUALOPSCTRIA_CONFIG_parallel_set");
    set_by_env(cfg_dualop.parallel_update,                "ESPRESO_DUALOPSCTRIA_CONFIG_parallel_update");
    set_by_env(cfg_dualop.parallel_apply,                 "ESPRESO_DUALOPSCTRIA_CONFIG_parallel_apply");
    set_by_env(cfg_dualop.mainloop_update_split,          "ESPRESO_DUALOPSCTRIA_CONFIG_mainloop_update_split");
    set_by_env(cfg_dualop.gpu_wait_after_mainloop_update, "ESPRESO_DUALOPSCTRIA_CONFIG_gpu_wait_after_mainloop_update");
    set_by_env(cfg_dualop.inner_timers,                   "ESPRESO_DUALOPSCTRIA_CONFIG_inner_timers");
    set_by_env(cfg_dualop.outer_timers,                   "ESPRESO_DUALOPSCTRIA_CONFIG_outer_timers");
    set_by_env(cfg_dualop.print_parameters,               "ESPRESO_DUALOPSCTRIA_CONFIG_print_parameters");

    replace_unset_configs<T,I>(cfg_sc, cfg_dualop);
}

template<typename T, typename I>
TotalFETIExplicitScTriaGpu<T,I>::TotalFETIExplicitScTriaGpu(FETI<T> &feti)
: DualOperator<T>(feti)
{
    if(!gpu::mgm::is_linked()) {
        eslog::error("TotalFETIExplicitScTriaGpu: not supported, espreso compiled without GPU support\n");
    }

    setup_configs<T,I>(op_sc_config, cfg);

    device = gpu::mgm::get_device_by_mpi(info::mpi::rank, info::mpi::size);
}



template<typename T, typename I>
TotalFETIExplicitScTriaGpu<T,I>::~TotalFETIExplicitScTriaGpu()
{
    gpu::mgm::memfree_device(ws_persistent);
    gpu::mgm::memfree_device(ws_tmp_for_cbmba);

    for(size_t i = 0; i < n_queues; i++) gpu::dnblas::handle_destroy(handles_dense[i]);
    for(size_t i = 0; i < n_queues; i++) gpu::spblas::handle_destroy(handles_sparse[i]);
    for(gpu::mgm::queue & q : queues) gpu::mgm::queue_destroy(q);
    gpu::mgm::queue_destroy(main_q);
}



template<typename T, typename I>
void TotalFETIExplicitScTriaGpu<T,I>::info()
{
    eslog::info(" = EXPLICIT TOTAL FETI OPERATOR USING TRIANGULAR SC ON GPU                                   = \n");
    eslog::info(" =   EXTERNAL SPARSE SOLVER               %50s = \n", DirectSparseSolver<T>::name());

    if(cfg.print_parameters) {
        auto order_to_string = [](char order){ switch(order){ case 'R': return "ROW_MAJOR"; case 'C': return "COL_MAJOR"; default: return "UNDEFINED"; }};
        auto bool_to_string = [](bool val){ return val ? "TRUE" : "FALSE";};
        auto loop_split_to_string = [](char val){ switch(val){ case 'C': return "COMBINED"; case 'S': return "SEPARATE"; default: return "UNDEFINED"; }};
        auto trsm_strategy_to_string = [](char val){ switch(val){ case 'R': return "SPLIT RHS TO BLOCK COLUMNS"; case 'F': return "SPLIT FACTOR TO BLOCKS"; default: return "UNDEFINED"; } };
        auto partalg_to_string = [](char alg){ switch(alg){ case 'U': return "UNIFORM"; case 'M': return "MINIMAL WORK"; default: return "UNDEFINED"; } };
        auto spdncrit_to_string = [](char val){ switch(val){ case 'S': return "SPARSE ONLY"; case 'D': return "DENSE ONLY"; case 'C': return "FRACTION OF CHUNKS IS SPARSE"; case 'Z': return "FRACTION ON SIZE IS SPARSE"; case 'T': return "DECIDE BASED ON DENSITY"; default: return "UNDEFINED"; } };
        auto spdn_to_string = [](char val){ switch(val){ case 'S': return "SPARSE"; case 'D': return "DENSE"; default: return "UNDEFINED"; } };
        auto prune_to_string = [](char val){ switch(val){ case 'N': return "NOTHING"; case 'R': return "PRUNE ROWS"; case 'C': return "PRUNE COLS"; case 'A': return "PRUNE_ALL"; default: return "UNDEFINED"; } };
        auto herk_strategy_to_string = [](char val){ switch(val){ case 'T': return "STAIRS"; case 'Q': return "SQUARES"; default: return "UNDEFINED"; } };

        eslog::info(" =   %-50s       %+30s = \n", "parallel_set", bool_to_string(cfg.parallel_set));
        eslog::info(" =   %-50s       %+30s = \n", "parallel_set", bool_to_string(cfg.parallel_update));
        eslog::info(" =   %-50s       %+30s = \n", "parallel_set", bool_to_string(cfg.parallel_apply));
        eslog::info(" =   %-50s       %+30s = \n", "mainloop_update_split", loop_split_to_string(cfg.mainloop_update_split));
        eslog::info(" =   %-50s       %+30s = \n", "gpu_wait_after_mainloop_update", bool_to_string(cfg.gpu_wait_after_mainloop_update));
        eslog::info(" =   %-50s       %+30s = \n", "order_F", order_to_string(cfg.order_F));
        eslog::info(" =   %-50s       %+30s = \n", "order_X", order_to_string(op_sc_config.order_X));
        eslog::info(" =   %-50s       %+30s = \n", "order_L", order_to_string(op_sc_config.order_L));
        eslog::info(" =   %-50s       %+30s = \n", "trsm.strategy", trsm_strategy_to_string(op_sc_config.cfg_trsm.strategy));
        eslog::info(" =   %-50s       %+30s = \n", "trsm.partition.algorithm", partalg_to_string(op_sc_config.cfg_trsm.partition.algorithm));
        eslog::info(" =   %-50s       %+30d = \n", "trsm.partition.parameter", op_sc_config.cfg_trsm.partition.parameter);
        eslog::info(" =   %-50s       %+30s = \n", "trsm.splitrhs.factor_order_sp", order_to_string(op_sc_config.cfg_trsm.splitrhs.factor_order_sp));
        eslog::info(" =   %-50s       %+30s = \n", "trsm.splitrhs.factor_order_dn", order_to_string(op_sc_config.cfg_trsm.splitrhs.factor_order_dn));
        eslog::info(" =   %-50s       %+30s = \n", "trsm.splitrhs.spdn_criteria", spdncrit_to_string(op_sc_config.cfg_trsm.splitrhs.spdn_criteria));
        eslog::info(" =   %-50s       %30.3f = \n", "trsm.splitrhs.spdn_param", op_sc_config.cfg_trsm.splitrhs.spdn_param);
        eslog::info(" =   %-50s       %+30s = \n", "trsm.splitfactor.trsm_factor_spdn", spdn_to_string(op_sc_config.cfg_trsm.splitfactor.trsm_factor_spdn));
        eslog::info(" =   %-50s       %+30s = \n", "trsm.splitfactor.trsm_factor_order", order_to_string(op_sc_config.cfg_trsm.splitfactor.trsm_factor_order));
        eslog::info(" =   %-50s       %+30s = \n", "trsm.splitfactor.gemm_factor_order_sp", order_to_string(op_sc_config.cfg_trsm.splitfactor.gemm_factor_order_sp));
        eslog::info(" =   %-50s       %+30s = \n", "trsm.splitfactor.gemm_factor_order_dn", order_to_string(op_sc_config.cfg_trsm.splitfactor.gemm_factor_order_dn));
        eslog::info(" =   %-50s       %+30s = \n", "trsm.splitfactor.gemm_factor_prune", prune_to_string(op_sc_config.cfg_trsm.splitfactor.gemm_factor_prune));
        eslog::info(" =   %-50s       %+30s = \n", "trsm.splitfactor.gemm_spdn_criteria", spdncrit_to_string(op_sc_config.cfg_trsm.splitfactor.gemm_spdn_criteria));
        eslog::info(" =   %-50s       %30.3f = \n", "trsm.splitfactor.gemm_spdn_param", op_sc_config.cfg_trsm.splitfactor.gemm_spdn_param);
        eslog::info(" =   %-50s       %+30s = \n", "herk.strategy", herk_strategy_to_string(op_sc_config.cfg_herk.strategy));
        eslog::info(" =   %-50s       %+30s = \n", "herk.partition_algorithm", partalg_to_string(op_sc_config.cfg_herk.partition_algorithm));
        eslog::info(" =   %-50s       %+30d = \n", "herk.partition_parameter", op_sc_config.cfg_herk.partition_parameter);
    }

    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_domain; }).to_string("  Domain volume [dofs]").c_str());
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_interface; }).to_string("  Domain surface [dofs]").c_str());
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}



template<typename T, typename I>
void TotalFETIExplicitScTriaGpu<T,I>::set(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitScTriaGpu::set");

    n_domains = feti.K.size();
    n_queues = omp_get_max_threads();
    
    gpu::mgm::init_gpu(device);
    gpu::mgm::set_device(device);

    queues.resize(n_queues);
    handles_dense.resize(n_queues);
    handles_sparse.resize(n_queues);

    domain_data.resize(n_domains);
    
    gpu::mgm::queue_create(main_q);
    for(gpu::mgm::queue & q : queues) gpu::mgm::queue_create(q);
    for(size_t i = 0; i < n_queues; i++) gpu::dnblas::handle_create(handles_dense[i], queues[i]);
    for(size_t i = 0; i < n_queues; i++) gpu::spblas::handle_create(handles_sparse[i], queues[i]);

    gpu::dnblas::init_library(main_q);
    gpu::mgm::queue_wait(main_q);

    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        data.n_dofs_domain = feti.B1[di].ncols;
        data.n_dofs_interface = feti.B1[di].nrows;
    }

    ssize_t free_mem_before_Falloc = gpu::mgm::get_device_memory_free();
    d_Fs_allocated.resize((n_domains - 1) / 2 + 1);
    {
        std::vector<size_t> domain_idxs_sorted_by_f_size_desc(n_domains);
        for(size_t di = 0; di < n_domains; di++) {
            domain_idxs_sorted_by_f_size_desc[di] = di;
        }
        std::sort(domain_idxs_sorted_by_f_size_desc.rbegin(), domain_idxs_sorted_by_f_size_desc.rend(), [&](size_t dl, size_t dr){ return domain_data[dl].n_dofs_interface < domain_data[dr].n_dofs_interface; });
        for(size_t i = 0; i < n_domains; i++) {
            size_t di = domain_idxs_sorted_by_f_size_desc[i];
            size_t di_bigger = domain_idxs_sorted_by_f_size_desc[(i / 2) * 2];
            size_t allocated_F_index = i / 2;
            MatrixDenseData_new<T> & d_F_allocd = d_Fs_allocated[allocated_F_index];
            per_domain_stuff & data_bigger = domain_data[di_bigger];
            per_domain_stuff & data_di = domain_data[di];
            if(di == di_bigger) {
                d_F_allocd.set(data_bigger.n_dofs_interface + 1, data_bigger.n_dofs_interface, cfg.order_F, AllocatorGPU_new::get_singleton());
                d_F_allocd.alloc();
                data_di.d_F = d_F_allocd.get_submatrix_view(1, data_di.n_dofs_interface + 1, 0, data_di.n_dofs_interface);
                data_di.d_F.prop.uplo = 'L';
            }
            else {
                data_di.d_F = d_F_allocd.get_submatrix_view(0, data_di.n_dofs_interface , 0, data_di.n_dofs_interface);
                data_di.d_F.prop.uplo = 'U';
            }
            data_di.d_F_old = MatrixDenseView_new<T>::template to_old<I,gpu::mgm::Ad>(data_di.d_F);
        }
    }
    ssize_t free_mem_after_Falloc = gpu::mgm::get_device_memory_free();
    ssize_t allocd_in_Falloc = free_mem_before_Falloc - free_mem_after_Falloc;
    stacktimer::info("TotalFETIExplicitScTriaGpu::set allocd_in_Falloc %zd", allocd_in_Falloc);

    stacktimer::push("TotalFETIExplicitScTriaGpu::set setup");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        gpu::mgm::queue & q = queues[di % n_queues];
        gpu::spblas::handle & hs = handles_sparse[di % n_queues];
        gpu::dnblas::handle & hd = handles_dense[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        stacktimer::info("TotalFETIExplicitScTriaGpu::set setup subdomain %zu", di);

        math::combine(data.Kreg, feti.K[di], feti.RegMat[di]);
        if constexpr(utils::is_real<T>())    data.Kreg.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
        if constexpr(utils::is_complex<T>()) data.Kreg.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE;
        data.Kreg.shape = feti.K[di].shape;

        data.solver_Kreg.commit(data.Kreg);
    
        stacktimer::push("set_symbolic_factorization");
        data.solver_Kreg.symbolicFactorization();
        data.n_nz_factor = data.solver_Kreg.getFactorNnz();
        stacktimer::pop();

        data.h_Bt = MatrixCsxView_new<T,I>::from_old(feti.B1[di]).get_transposed_reordered_view();

        data.op_sc = std::make_unique<gpu::operations::sc_symm_hcsx_ddny_tria<T,I>>();
        data.op_sc->set_config(op_sc_config);
        data.op_sc->set_handles(q, hs, hd);
        data.op_sc->set_coefficients(-1);
        data.op_sc->set_h_A11_solver(&data.solver_Kreg);
        data.op_sc->set_h_A12(&data.h_Bt);
        data.op_sc->set_d_sc(&data.d_F);
        data.op_sc->setup();

        data.d_apply_x.resize(data.n_dofs_interface);
        data.d_apply_y.resize(data.n_dofs_interface);
        data.d_apply_z.resize(data.n_dofs_domain);
        data.d_apply_w.resize(data.n_dofs_domain);
        data.d_applyg_D2C.resize(data.n_dofs_interface);
    }
    if(!cfg.inner_timers) stacktimer::enable();
    stacktimer::pop();

    // clean up the mess from buggy openmp in clang
    utils::run_dummy_parallel_region();

    total_wss_internal = 0;
    total_wss_persistent = 0;
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        total_wss_internal += data.op_sc->get_wss_internal();
        total_wss_persistent += data.op_sc->get_wss_persistent();
    }
    total_wss_internal = utils::round_up((total_wss_internal * 105) / 100, gpu::mgm::get_natural_pitch_align());

    stacktimer::info("TotalFETIExplicitScTriaGpu::set total_wss_internal %zu", total_wss_internal);
    stacktimer::info("TotalFETIExplicitScTriaGpu::set total_wss_persistent %zu", total_wss_persistent);
    stacktimer::info("TotalFETIExplicitScTriaGpu::set max_wss_tmp_preprocess %zu", std::max_element(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & l, const per_domain_stuff & r){ return l.op_sc->get_wss_tmp_preprocess() < r.op_sc->get_wss_tmp_preprocess(); })->op_sc->get_wss_tmp_preprocess());
    stacktimer::info("TotalFETIExplicitScTriaGpu::set max_wss_tmp_perform %zu", std::max_element(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & l, const per_domain_stuff & r){ return l.op_sc->get_wss_tmp_perform() < r.op_sc->get_wss_tmp_perform(); })->op_sc->get_wss_tmp_perform());

    ws_persistent = gpu::mgm::memalloc_device(total_wss_persistent);
    ator_ws_persistent = std::make_unique<AllocatorArena_new>(false, true, gpu::mgm::get_natural_pitch_align());
    ator_ws_persistent->set(ws_persistent, total_wss_persistent);
    
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        data.op_sc->set_ws_persistent(ator_ws_persistent->alloc(data.op_sc->get_wss_persistent()));
    }

    {
        d_applyg_x_cluster.resize(feti.lambdas.size);
        d_applyg_y_cluster.resize(feti.lambdas.size);
        Vector_Dense<T*,I,gpu::mgm::Ah> h_applyg_xs_pointers;
        Vector_Dense<T*,I,gpu::mgm::Ah> h_applyg_ys_pointers;
        Vector_Dense<I,I,gpu::mgm::Ah> h_applyg_n_dofs_interfaces;
        Vector_Dense<I*,I,gpu::mgm::Ah> h_applyg_D2Cs_pointers;
        h_applyg_xs_pointers.resize(n_domains);
        h_applyg_ys_pointers.resize(n_domains);
        h_applyg_n_dofs_interfaces.resize(n_domains);
        h_applyg_D2Cs_pointers.resize(n_domains);
        d_applyg_xs_pointers.resize(n_domains);
        d_applyg_ys_pointers.resize(n_domains);
        d_applyg_n_dofs_interfaces.resize(n_domains);
        d_applyg_D2Cs_pointers.resize(n_domains);
        for(size_t di = 0; di < n_domains; di++) {
            h_applyg_xs_pointers.vals[di] = domain_data[di].d_apply_x.vals;
            h_applyg_ys_pointers.vals[di] = domain_data[di].d_apply_y.vals;
            h_applyg_n_dofs_interfaces.vals[di] = domain_data[di].n_dofs_interface;
            h_applyg_D2Cs_pointers.vals[di] = domain_data[di].d_applyg_D2C.vals;
        }
        gpu::mgm::copy_submit(main_q, d_applyg_xs_pointers,       h_applyg_xs_pointers);
        gpu::mgm::copy_submit(main_q, d_applyg_ys_pointers,       h_applyg_ys_pointers);
        gpu::mgm::copy_submit(main_q, d_applyg_n_dofs_interfaces, h_applyg_n_dofs_interfaces);
        gpu::mgm::copy_submit(main_q, d_applyg_D2Cs_pointers,     h_applyg_D2Cs_pointers);
        for(size_t di = 0; di < n_domains; di++) {
            gpu::mgm::copy_submit(main_q, domain_data[di].d_applyg_D2C.vals, feti.D2C[di].data(), feti.D2C[di].size());
        }
    }

    size_t free_mem = gpu::mgm::get_device_memory_free();
    size_t mem_capacity = gpu::mgm::get_device_memory_capacity();
    size_t reserve = (mem_capacity * 5) / 100;
    wss_tmp_for_cbmba = utils::round_down(free_mem - reserve - total_wss_internal, gpu::mgm::get_natural_pitch_align());
    ws_tmp_for_cbmba = gpu::mgm::memalloc_device(wss_tmp_for_cbmba);
    ator_tmp_cbmba = std::make_unique<AllocatorCBMB_new>(false, true, gpu::mgm::get_natural_pitch_align(), ws_tmp_for_cbmba, wss_tmp_for_cbmba);

    stacktimer::info("TotalFETIExplicitScTriaGpu::set cbmba_capacity %zu", wss_tmp_for_cbmba);

    ssize_t free_mem_before_preprocess = gpu::mgm::get_device_memory_free();

    stacktimer::push("TotalFETIExplicitScTriaGpu::set preprocess");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        gpu::mgm::queue & q = queues[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        stacktimer::info("TotalFETIExplicitScTriaGpu::set preprocess subdomain %zu", di);

        void * ws_tmp = ator_tmp_cbmba->alloc(data.op_sc->get_wss_tmp_preprocess());

        data.op_sc->preprocess_submit(ws_tmp);

        gpu::mgm::submit_host_function(q, [&,ws_tmp,di](){
            void * ws_tmp_ = ws_tmp;
            ator_tmp_cbmba->free(ws_tmp_);
        });
    }
    gpu::mgm::device_wait();
    if(!cfg.inner_timers) stacktimer::enable();
    stacktimer::pop();

    ssize_t free_mem_after_preprocess = gpu::mgm::get_device_memory_free();
    ssize_t allocd_in_preprocess = free_mem_before_preprocess - free_mem_after_preprocess;
    stacktimer::info("TotalFETIExplicitScTriaGpu::set allocd_in_preprocess %zd", allocd_in_preprocess);

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIExplicitScTriaGpu<T,I>::update(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitScTriaGpu::update");

    gpu::mgm::set_device(device);

    auto loop_part_1_factorize = [this](size_t di) {
        per_domain_stuff & data = domain_data[di];

        stacktimer::info("TotalFETIExplicitScTriaGpu::update subdomain %zu", di);

        math::sumCombined(data.Kreg, T{1.0}, feti.K[di], feti.RegMat[di]);

        stacktimer::push("update_commit_matrix_to_solver");
        data.solver_Kreg.commit(data.Kreg);
        stacktimer::pop();

        stacktimer::push("update_numerical_factorization");
        data.solver_Kreg.numericalFactorization();
        stacktimer::pop();
    };
    auto loop_part_2_assemble = [this](size_t di) {
        gpu::mgm::queue & q = queues[di % n_queues];
        per_domain_stuff & data = domain_data[di];
        
        stacktimer::push("update_cbmba_allocation");
        void * ws_tmp = ator_tmp_cbmba->alloc(data.op_sc->get_wss_tmp_perform());
        stacktimer::pop();

        data.op_sc->perform_submit(ws_tmp);

        gpu::mgm::submit_host_function(q, [&,ws_tmp](){
            void * ws_tmp_ = ws_tmp;
            ator_tmp_cbmba->free(ws_tmp_);
        });
    };

    stacktimer::push("update_mainloop");
    if(cfg.mainloop_update_split == 'C') {
        stacktimer::push("update_mainloop_combined");
        if(!cfg.inner_timers) stacktimer::disable();
        #pragma omp parallel for schedule(static,1) if(cfg.parallel_update)
        for(size_t di = 0; di < n_domains; di++) {
            loop_part_1_factorize(di);
            loop_part_2_assemble(di);
        };
        if(!cfg.inner_timers) stacktimer::enable();
        stacktimer::pop();

        // clean up the mess from buggy openmp in clang
        utils::run_dummy_parallel_region();

        if(cfg.gpu_wait_after_mainloop_update) {
            stacktimer::push("update_wait_to_finish_after_mainloop");
            gpu::mgm::device_wait();
            stacktimer::pop();
        }
    }
    if(cfg.mainloop_update_split == 'S') {
        stacktimer::push("update_mainloop_sepatare_factorize");
        if(!cfg.inner_timers) stacktimer::disable();
        #pragma omp parallel for schedule(static,1) if(cfg.parallel_update)
        for(size_t di = 0; di < n_domains; di++) {
            loop_part_1_factorize(di);
        }
        if(!cfg.inner_timers) stacktimer::enable();
        stacktimer::pop();

        // clean up the mess from buggy openmp in clang
        utils::run_dummy_parallel_region();

        stacktimer::push("update_mainloop_sepatare_assemble");
        stacktimer::push("update_mainloop_sepatare_assemble_submit");
        if(!cfg.inner_timers) stacktimer::disable();
        #pragma omp parallel for schedule(static,1) if(cfg.parallel_update)
        for(size_t di = 0; di < n_domains; di++) {
            loop_part_2_assemble(di);
        }
        if(!cfg.inner_timers) stacktimer::enable();
        stacktimer::pop();

        if(cfg.gpu_wait_after_mainloop_update) {
            stacktimer::push("update_mainloop_sepatare_assemble_wait");
            gpu::mgm::device_wait();
            stacktimer::pop();
        }
        stacktimer::pop();
    }
    stacktimer::pop();

    // clean up the mess from buggy openmp in clang
    utils::run_dummy_parallel_region();

    stacktimer::push("update_compute_vector_d");
    {
        if (feti.updated.B) {
            d.resize();
        }
        // just use the cpu solver
        std::vector<Vector_Dense<T,I>> Kplus_fs(n_domains);
        #pragma omp parallel for schedule(static,1) if(cfg.parallel_update)
        for(size_t di = 0; di < n_domains; di++) {
            domain_data[di].solver_Kreg.solve(feti.f[di], Kplus_fs[di]);
        }
        applyB(feti, Kplus_fs, d);
        d.synchronize();
        math::add(d, T{-1}, feti.c);
    }
    stacktimer::pop();

    stacktimer::push("update_final_wait");
    gpu::mgm::device_wait();
    stacktimer::pop();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIExplicitScTriaGpu<T,I>::_apply(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitScTriaGpu::apply");

    gpu::mgm::set_device(device);

    // copy x_cluster to device
    gpu::mgm::copy_submit(main_q, d_applyg_x_cluster, x_cluster);

    // scatter
    gpu::kernels::DCmap_scatter(main_q, d_applyg_xs_pointers, d_applyg_n_dofs_interfaces, d_applyg_x_cluster, d_applyg_D2Cs_pointers);

    gpu::mgm::queue_async_barrier({main_q}, queues);

    // apply
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_apply)
    for(size_t di = 0; di < n_domains; di++) {
        gpu::dnblas::handle & hd = handles_dense[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        gpu::dnblas::hemv<T,I>(hd, data.d_F.nrows, data.d_F.vals, data.d_F.ld, data.d_F.order, 'N', data.d_F.prop.uplo, data.d_apply_x.vals, data.d_apply_y.vals);
    }

    // zerofill y_cluster on device
    gpu::mgm::memset_submit(main_q, d_applyg_y_cluster.vals, d_applyg_y_cluster.size * sizeof(T), 0);

    gpu::mgm::queue_async_barrier(queues, {main_q});

    // gather
    gpu::kernels::DCmap_gather(main_q, d_applyg_ys_pointers, d_applyg_n_dofs_interfaces, d_applyg_y_cluster, d_applyg_D2Cs_pointers);

    // copy y_cluster from device
    gpu::mgm::copy_submit(main_q, y_cluster, d_applyg_y_cluster);

    // wait
    gpu::mgm::device_wait();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template <typename T, typename I>
void TotalFETIExplicitScTriaGpu<T,I>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    _apply(x, y);
    y.synchronize();
}



template <typename T, typename I>
void TotalFETIExplicitScTriaGpu<T,I>::apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y)
{
    Vector_Dual<T> _x, _y;
    _x.size = _y.size = x.ncols;
    for (int r = 0; r < x.nrows; ++r) {
        _x.vals = x.vals + x.ncols * r;
        _y.vals = y.vals + y.ncols * r;
        _apply(_x, _y);
    }
    y.synchronize();
}



template<typename T, typename I>
void TotalFETIExplicitScTriaGpu<T,I>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    #pragma omp parallel for schedule(static,1)
    for (size_t di = 0; di < n_domains; ++di) {
        Vector_Dense<T,I> z;
        z.resize(y[di]);
        applyBt(feti, di, x, z, T{-1});
        math::add(z, T{1}, feti.f[di]);
        domain_data[di].solver_Kreg.solve(z, y[di]);
    }
}



template<typename T, typename I>
void TotalFETIExplicitScTriaGpu<T,I>::print(const step::Step &step)
{
    eslog::error("TotalFETIExplicitScTriaGpu::print not implemented");
}





#define INSTANTIATE_T_I(T,I) \
template class TotalFETIExplicitScTriaGpu<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        INSTANTIATE_T(std::complex<double>)

    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I

}
