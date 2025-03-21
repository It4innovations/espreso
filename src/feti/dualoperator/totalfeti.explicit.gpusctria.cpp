
#include "totalfeti.explicit.gpusctria.h"
#include "math/wrappers/math.blas.h"
#include "feti/common/applyB.h"
#include "basis/utilities/minmaxavg.h"
#include "my_timer.h"
#include "gpu/gpu_kernels.h"

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

template<typename T, typename I>
static void setup_sc_config(typename gpu::operations::sc_symm_hcsx_ddny_tria<T,I>::config & cfg)
{
    set_by_env(cfg.order_X,                                   "ESPRESO_DUALOPGPUSCTRIA_CONFIG_order_X");
    set_by_env(cfg.order_L,                                   "ESPRESO_DUALOPGPUSCTRIA_CONFIG_order_L");
    set_by_env(cfg.cfg_trsm.strategy,                         "ESPRESO_DUALOPGPUSCTRIA_CONFIG_trsm_strategy");
    set_by_env(cfg.cfg_trsm.partition.algorithm,              "ESPRESO_DUALOPGPUSCTRIA_CONFIG_trsm_partition_algorithm");
    set_by_env(cfg.cfg_trsm.partition.parameter,              "ESPRESO_DUALOPGPUSCTRIA_CONFIG_trsm_partition_parameter");
    set_by_env(cfg.cfg_trsm.splitrhs.factor_order_sp,         "ESPRESO_DUALOPGPUSCTRIA_CONFIG_trsm_splitrhs_factor_order_sp");
    set_by_env(cfg.cfg_trsm.splitrhs.factor_order_dn,         "ESPRESO_DUALOPGPUSCTRIA_CONFIG_trsm_splitrhs_factor_order_dn");
    set_by_env(cfg.cfg_trsm.splitrhs.spdn_criteria,           "ESPRESO_DUALOPGPUSCTRIA_CONFIG_trsm_splitrhs_spdn_criteria");
    set_by_env(cfg.cfg_trsm.splitrhs.spdn_param,              "ESPRESO_DUALOPGPUSCTRIA_CONFIG_trsm_splitrhs_spdn_param");
    set_by_env(cfg.cfg_trsm.splitfactor.trsm_factor_spdn,     "ESPRESO_DUALOPGPUSCTRIA_CONFIG_trsm_splitfactor_trsm_factor_spdn");
    set_by_env(cfg.cfg_trsm.splitfactor.trsm_factor_order,    "ESPRESO_DUALOPGPUSCTRIA_CONFIG_trsm_splitfactor_trsm_factor_order");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_factor_order_sp, "ESPRESO_DUALOPGPUSCTRIA_CONFIG_trsm_splitfactor_gemm_factor_order_sp");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_factor_order_dn, "ESPRESO_DUALOPGPUSCTRIA_CONFIG_trsm_splitfactor_gemm_factor_order_dn");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_factor_prune,    "ESPRESO_DUALOPGPUSCTRIA_CONFIG_trsm_splitfactor_gemm_factor_prune");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_spdn_criteria,   "ESPRESO_DUALOPGPUSCTRIA_CONFIG_trsm_splitfactor_gemm_spdn_criteria");
    set_by_env(cfg.cfg_trsm.splitfactor.gemm_spdn_param,      "ESPRESO_DUALOPGPUSCTRIA_CONFIG_trsm_splitfactor_gemm_spdn_param");
    set_by_env(cfg.cfg_herk.strategy,                         "ESPRESO_DUALOPGPUSCTRIA_CONFIG_herk_strategy");
    set_by_env(cfg.cfg_herk.partition_algorithm,              "ESPRESO_DUALOPGPUSCTRIA_CONFIG_herk_partition_algorithm");
    set_by_env(cfg.cfg_herk.partition_parameter,              "ESPRESO_DUALOPGPUSCTRIA_CONFIG_herk_partition_parameter");
}

template<typename T, typename I>
TotalFETIExplicitGpuScTria<T,I>::TotalFETIExplicitGpuScTria(FETI<T> &feti)
: DualOperator<T>(feti)
{
    if(!gpu::mgm::is_linked()) {
        eslog::error("TotalFETIExplicitGpuScTria: not supported, espreso compiled without GPU support\n");
    }
    
    device = gpu::mgm::get_device_by_mpi(info::mpi::rank, info::mpi::size);
}



template<typename T, typename I>
TotalFETIExplicitGpuScTria<T,I>::~TotalFETIExplicitGpuScTria()
{
    gpu::mgm::memfree_device(ws_persistent);
    gpu::mgm::memfree_device(ws_tmp_for_cbmba);

    for(size_t i = 0; i < n_queues; i++) gpu::dnblas::handle_destroy(handles_dense[i]);
    for(size_t i = 0; i < n_queues; i++) gpu::spblas::handle_destroy(handles_sparse[i]);
    for(gpu::mgm::queue & q : queues) gpu::mgm::queue_destroy(q);
    gpu::mgm::queue_destroy(main_q);
}



template<typename T, typename I>
void TotalFETIExplicitGpuScTria<T,I>::info()
{
    eslog::info(" = EXPLICIT TOTAL FETI OPERATOR ON GPU USING TRIANGULAR SC                                   = \n");
    eslog::info(" =   EXTERNAL SPARSE SOLVER               %50s = \n", DirectSparseSolver<T>::name());
    // eslog::info(" =   %-50s       %+30s = \n", "APPLY_WHERE", apply_on_gpu ? "GPU" : "CPU");
    // eslog::info(minmaxavg<double>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.F.nrows * data.F.get_ld() * sizeof(T) / (1024.0 * 1024.0); }).to_string("  F MEMORY [MB]").c_str());
    // eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_domain; }).to_string("  Domain volume [dofs]").c_str());
    // eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_interface; }).to_string("  Domain surface [dofs]").c_str());
    // eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}



template<typename T, typename I>
void TotalFETIExplicitGpuScTria<T,I>::set(const step::Step &step)
{
    char order_F = '_';
    typename gpu::operations::sc_symm_hcsx_ddny_tria<T,I>::config op_sc_config;
    setup_sc_config<T,I>(op_sc_config);
    set_by_env(order_F, "ESPRESO_DUALOPGPUSCTRIA_CONFIG_order_F");
    set_by_env(parallel_set, "ESPRESO_DUALOPGPUSCTRIA_CONFIG_parallel_set");
    set_by_env(parallel_update, "ESPRESO_DUALOPGPUSCTRIA_CONFIG_parallel_update");
    set_by_env(parallel_apply, "ESPRESO_DUALOPGPUSCTRIA_CONFIG_parallel_apply");

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
                d_F_allocd.set(data_bigger.n_dofs_interface + 1, data_bigger.n_dofs_interface, order_F, AllocatorGPU_new::get_singleton());
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

    #pragma omp parallel for schedule(static,1) if(parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        gpu::mgm::queue & q = queues[di % n_queues];
        gpu::spblas::handle & hs = handles_sparse[di % n_queues];
        gpu::dnblas::handle & hd = handles_dense[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        // Kreg = K + RegMat symbolic pattern
        // tm_Kreg_combine.start();
        {
            math::combine(data.Kreg, feti.K[di], feti.RegMat[di]);
            if constexpr(utils::is_real<T>())    data.Kreg.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
            if constexpr(utils::is_complex<T>()) data.Kreg.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE;
            data.Kreg.shape = feti.K[di].shape;
        }
        // tm_Kreg_combine.stop();

        // commit Kreg to solver (just symbolic pattern present now)
        // tm_solver_commit.start();
        {
            data.solver_Kreg.commit(data.Kreg);
        }
        // tm_solver_commit.stop();

        // symbolic factorization
        // tm_fact_symbolic.start();
        {
            data.solver_Kreg.symbolicFactorization();
            data.n_nz_factor = data.solver_Kreg.getFactorNnz();
        }
        // tm_fact_symbolic.stop();

        data.h_Bt = MatrixCsxView_new<T,I>::from_old(feti.B0[di]).get_transposed_reordered_view();

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

    total_wss_internal = 0;
    total_wss_persistent = 0;
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        total_wss_internal += data.op_sc->get_wss_internal();
        total_wss_persistent += data.op_sc->get_wss_persistent();
    }
    total_wss_internal = utils::round_up((total_wss_internal * 105) / 100, gpu::mgm::get_natural_pitch_align());

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

    #pragma omp parallel for schedule(static,1) if(parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        gpu::mgm::queue & q = queues[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        void * ws_tmp = ator_tmp_cbmba->alloc(data.op_sc->get_wss_tmp_preprocess());

        data.op_sc->preprocess_submit(ws_tmp);

        gpu::mgm::submit_host_function(q, [&,ws_tmp](){
            void * ws_tmp_ = ws_tmp;
            ator_tmp_cbmba->free(ws_tmp_);
        });
    }

    gpu::mgm::device_wait();
}



template<typename T, typename I>
void TotalFETIExplicitGpuScTria<T,I>::update(const step::Step &step)
{
    gpu::mgm::set_device(device);
    
    #pragma omp parallel for schedule(static,1) if(parallel_update)
    for(size_t di = 0; di < n_domains; di++) {
        gpu::mgm::queue & q = queues[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        // Kreg = K + RegMat numeric values
        // tm_Kreg_combine.start();
        {
            math::sumCombined(data.Kreg, T{1.0}, feti.K[di], feti.RegMat[di]);
        }
        // tm_Kreg_combine.stop();

        // commit Kreg to solver (with numeric values)
        // tm_solver_commit.start();
        {
            data.solver_Kreg.commit(data.Kreg);
        }
        // tm_solver_commit.stop();
        
        // numeric factorization
        // tm_fact_numeric.start();
        {
            data.solver_Kreg.numericalFactorization();
        }
        // tm_fact_numeric.stop();

        void * ws_tmp = ator_tmp_cbmba->alloc(data.op_sc->get_wss_tmp_preprocess());

        data.op_sc->perform_submit(ws_tmp);

        gpu::mgm::submit_host_function(q, [&,ws_tmp](){
            void * ws_tmp_ = ws_tmp;
            ator_tmp_cbmba->free(ws_tmp_);
        });
    }

    gpu::mgm::device_wait();
}



template<typename T, typename I>
void TotalFETIExplicitGpuScTria<T,I>::_apply(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    gpu::mgm::set_device(device);

    // copy x_cluster to device
    // tm_copyin.start();
    gpu::mgm::copy_submit(main_q, d_applyg_x_cluster, x_cluster);
    // if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
    // tm_copyin.stop();

    // scatter
    // tm_scatter.start();
    gpu::kernels::DCmap_scatter(main_q, d_applyg_xs_pointers, d_applyg_n_dofs_interfaces, d_applyg_x_cluster, d_applyg_D2Cs_pointers);
    // if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
    // tm_scatter.stop();

    gpu::mgm::queue_async_barrier({main_q}, queues);

    // apply
    // tm_mv_outer.start();
    #pragma omp parallel for schedule(static,1) if(parallel_apply)
    for(size_t di = 0; di < n_domains; di++) {
        // gpu::mgm::queue & q = queues[di % n_queues];
        gpu::dnblas::handle & hd = handles_dense[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        // tm_mv.start();
        gpu::dnblas::hemv<T,I>(hd, data.d_F.nrows, data.d_F.vals, data.d_F.ld, data.d_F.order, 'N', data.d_F.prop.uplo, data.d_apply_x.vals, data.d_apply_y.vals);
        // if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        // tm_mv.stop();
    }
    // tm_mv_outer.stop();

    // zerofill y_cluster on device
    // tm_zerofill.start();
    gpu::mgm::memset_submit(main_q, d_applyg_y_cluster.vals, d_applyg_y_cluster.size * sizeof(T), 0);
    // if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
    // tm_zerofill.stop();

    gpu::mgm::queue_async_barrier(queues, {main_q});

    // gather
    // tm_gather.start();
    gpu::kernels::DCmap_gather(main_q, d_applyg_ys_pointers, d_applyg_n_dofs_interfaces, d_applyg_y_cluster, d_applyg_D2Cs_pointers);
    // if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
    // tm_gather.stop();

    // copy y_cluster from device
    // tm_copyout.start();
    gpu::mgm::copy_submit(main_q, y_cluster, d_applyg_y_cluster);
    // if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
    // tm_copyout.stop();

    // wait
    // tm_wait.start();
    gpu::mgm::device_wait();
    // tm_wait.stop();

    // tm_total.stop();
}



template <typename T, typename I>
void TotalFETIExplicitGpuScTria<T,I>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    _apply(x, y);
    y.synchronize();
}



template <typename T, typename I>
void TotalFETIExplicitGpuScTria<T,I>::apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y)
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
void TotalFETIExplicitGpuScTria<T,I>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
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
void TotalFETIExplicitGpuScTria<T,I>::print(const step::Step &step)
{
    eslog::error("TotalFETIExplicitGpuScTria::print not implemented");
}





#define INSTANTIATE_T_I(T,I) \
template class TotalFETIExplicitGpuScTria<T,I>;

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
