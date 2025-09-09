
#include "feti/dualoperator/totalfeti.explicit.generalschur.gpu.h"

#include "feti/common/applyB.h"
#include "basis/utilities/minmaxavg.h"
#include "basis/utilities/stacktimer.h"
#include "esinfo/meshinfo.h"
#include "gpu/gpu_kernels.h"

#include <algorithm>



namespace espreso {



template<typename T, typename I>
TotalFETIExplicitGeneralSchurGpu<T,I>::TotalFETIExplicitGeneralSchurGpu(FETI<T> &feti) : DualOperator<T>(feti)
    , main_q(feti.main_q)
    , queues(feti.queues)
    , handles_dense(feti.handles_dense)
    , handles_sparse(feti.handles_sparse)
{
    setup_config(cfg, feti.configuration);
}



template<typename T, typename I>
TotalFETIExplicitGeneralSchurGpu<T,I>::~TotalFETIExplicitGeneralSchurGpu()
{
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurGpu<T,I>::info()
{
    eslog::info(" = EXPLICIT TOTAL FETI OPERATOR USING GENERAL SCHUR ON GPU                                   = \n");

    if(cfg.print_config) {
        auto order_to_string = [](char order){ switch(order){ case 'R': return "ROW_MAJOR"; case 'C': return "COL_MAJOR"; default: return "UNDEFINED"; }};
        auto bool_to_string = [](bool val){ return val ? "TRUE" : "FALSE";};
        auto loop_split_to_string = [](char val){ switch(val){ case 'C': return "COMBINED"; case 'S': return "SEPARATE"; default: return "UNDEFINED"; }};
        auto schur_impl_to_string = [](schur_impl_t schur_impl){ switch(schur_impl) { case schur_impl_t::autoselect: return "autoselect"; case schur_impl_t::manual_simple: return "manual_simple"; case schur_impl_t::triangular: return "triangular"; default: return "UNDEFINED"; } };
        auto schur_impl_to_string_actual = [](schur_impl_t schur_impl){ return gpu::operations::schur_hcsx_ddny<T,I>::make(schur_impl)->get_name(); };

        eslog::info(" =   %-50s       %+30s = \n", "parallel_set", bool_to_string(cfg.parallel_set));
        eslog::info(" =   %-50s       %+30s = \n", "parallel_set", bool_to_string(cfg.parallel_update));
        eslog::info(" =   %-50s       %+30s = \n", "parallel_set", bool_to_string(cfg.parallel_apply));
        eslog::info(" =   %-50s       %+30s = \n", "mainloop_update_split", loop_split_to_string(cfg.mainloop_update_split));
        eslog::info(" =   %-50s       %+30s = \n", "synchronize_after_update_mainloops", loop_split_to_string(cfg.gpu_wait_after_mainloop_update));
        eslog::info(" =   %-50s       %+30s = \n", "inner_timers", bool_to_string(cfg.inner_timers));
        eslog::info(" =   %-50s       %+30s = \n", "outer_timers", bool_to_string(cfg.outer_timers));
        eslog::info(" =   %-50s       %+30s = \n", "order_F", order_to_string(cfg.order_F));
        eslog::info(" =   %-50s       %+30s = \n", "schur implementation", schur_impl_to_string(cfg.schur_impl));
        eslog::info(" =   %-50s       %+30s = \n", "schur implementation actual", schur_impl_to_string_actual(cfg.schur_impl));
    }

    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_domain; }).to_string("  Domain volume [dofs]").c_str());
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_interface; }).to_string("  Domain interface [dofs]").c_str());
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurGpu<T,I>::setup()
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralSchurGpu::setup");

    total_wss_gpu_persistent = 0;
    total_wss_gpu_internal = 0;

    ator_ws_gpu_persistent = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());

    n_domains = feti.K.size();
    n_queues = queues.size();

    domain_data.resize(n_domains);
    
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
                d_F_allocd.set(data_bigger.n_dofs_interface + (cfg.order_F == 'R'), data_bigger.n_dofs_interface + (cfg.order_F == 'C'), cfg.order_F, ator_ws_gpu_persistent.get());
            }
            data_di.d_F.set_view(data_di.n_dofs_interface, data_di.n_dofs_interface, d_F_allocd.ld, d_F_allocd.order, nullptr, d_F_allocd.ator);
            data_di.op_sub_F_from_allocd.set_matrix_src(&d_F_allocd);
            data_di.op_sub_F_from_allocd.set_matrix_dst(&data_di.d_F);
            if(di == di_bigger) {
                data_di.op_sub_F_from_allocd.set_bounds(0, data_di.d_F.nrows, 0, data_di.d_F.ncols);
                data_di.d_F.prop.uplo = (d_F_allocd.order == 'R') ? 'U' : 'L';
            }
            else {
                data_di.op_sub_F_from_allocd.set_bounds((int)(cfg.order_F == 'R'), data_di.d_F.nrows + (int)(cfg.order_F == 'R'), (int)(cfg.order_F == 'C'), data_di.d_F.ncols + (int)(cfg.order_F == 'C'));
                data_di.d_F.prop.uplo = (d_F_allocd.order == 'R') ? 'L' : 'U';
            }
            data_di.d_F.prop.symm = MatrixSymmetry_new::hermitian;
        }
    }

    size_t wss_pers_Fs = 0;
    for(const auto & d_F : d_Fs_allocated) {
        wss_pers_Fs += d_F.get_memory_impact();
    }
    total_wss_gpu_persistent += wss_pers_Fs;

    {
        std::vector<MatrixDenseView_new<T>*> Fs_vector(n_domains);
        for(size_t di = 0; di < n_domains; di++) {
            Fs_vector[di] = &domain_data[di].d_F;
        }

        applicator.set_config(cfg.parallel_apply, cfg.inner_timers);
        applicator.set_handles(&main_q, &queues, &handles_dense);
        applicator.set_dimensions(feti);
        applicator.set_memory('C', 'G');
        applicator.set_D2C_map(&feti.D2C);
        applicator.set_Fs(Fs_vector);
        applicator.set_apply_target(cfg.apply_where);
        applicator.setup();
    }
    total_wss_gpu_persistent += applicator.get_wss_gpu_persistent();

    stacktimer::push("TotalFETIExplicitGeneralSchurGpu::setup op_sc");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        gpu::mgm::queue & q = queues[di % n_queues];
        gpu::spblas::handle & hs = handles_sparse[di % n_queues];
        gpu::dnblas::handle & hd = handles_dense[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        stacktimer::info("TotalFETIExplicitGeneralSchurGpu::setup op_sc subdomain %zu", di);

        math::combine(data.Kreg_old, feti.K[di], feti.RegMat[di]);
        if constexpr(utils::is_real<T>())    data.Kreg_old.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
        if constexpr(utils::is_complex<T>()) data.Kreg_old.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE;
        data.Kreg_old.shape = feti.K[di].shape;

        data.Bt = MatrixCsxView_new<T,I>::from_old(feti.B1[di]).get_transposed_reordered_view();
        data.Kreg = MatrixCsxView_new<T,I>::from_old(data.Kreg_old);

        data.op_sc = gpu::operations::schur_hcsx_ddny<T,I>::make(cfg.schur_impl);
        data.op_sc->set_handles(q, hs, hd);
        data.op_sc->set_coefficients(-1);
        data.op_sc->set_matrix(&data.Kreg, &data.Bt, nullptr, nullptr);
        data.op_sc->set_sc(&data.d_F);
        data.op_sc->set_need_solve_A11(true);
        data.op_sc->setup();
    }
    if(!cfg.inner_timers) stacktimer::enable();
    stacktimer::pop();

    // clean up the mess from buggy openmp in clang
    utils::run_dummy_parallel_region();

    size_t op_sc_wss_gpu_internal = 0;
    size_t op_sc_wss_gpu_persistent = 0;
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        op_sc_wss_gpu_internal += data.op_sc->get_wss_internal();
        op_sc_wss_gpu_persistent += data.op_sc->get_wss_persistent();
    }
    op_sc_wss_gpu_internal += utils::round_up((op_sc_wss_gpu_internal * 105) / 100, gpu::mgm::get_natural_pitch_align());
    total_wss_gpu_internal += op_sc_wss_gpu_internal;
    total_wss_gpu_persistent += op_sc_wss_gpu_persistent;

    stacktimer::info("TotalFETIExplicitGeneralSchurGpu::setup wss_gpu_persistent_Fs %zd", wss_pers_Fs);
    stacktimer::info("TotalFETIExplicitGeneralSchurGpu::setup wss_gpu_persistent_applicator %zd", applicator.get_wss_gpu_persistent());
    stacktimer::info("TotalFETIExplicitGeneralSchurGpu::setup wss_gpu_persistent_op_sc %zu", op_sc_wss_gpu_persistent);
    stacktimer::info("TotalFETIExplicitGeneralSchurGpu::setup wss_gpu_internal_op_sc %zu", op_sc_wss_gpu_internal);
    stacktimer::info("TotalFETIExplicitGeneralSchurGpu::setup max_wss_tmp_preprocess_op_sc %zu", std::max_element(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & l, const per_domain_stuff & r){ return l.op_sc->get_wss_tmp_preprocess() < r.op_sc->get_wss_tmp_preprocess(); })->op_sc->get_wss_tmp_preprocess());
    stacktimer::info("TotalFETIExplicitGeneralSchurGpu::setup max_wss_tmp_perform_op_sc %zu", std::max_element(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & l, const per_domain_stuff & r){ return l.op_sc->get_wss_tmp_perform() < r.op_sc->get_wss_tmp_perform(); })->op_sc->get_wss_tmp_perform());

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurGpu<T,I>::set(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralSchurGpu::set");

    if(ws_gpu_persistent == nullptr) eslog::error("persistent gpu workspace is not set\n");

    ator_ws_gpu_persistent->set(ws_gpu_persistent, total_wss_gpu_persistent);

    for(auto & d_F : d_Fs_allocated) d_F.alloc();

    if(!cfg.inner_timers) stacktimer::disable();
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & my = domain_data[di];
        my.op_sub_F_from_allocd.perform();
    }
    if(!cfg.inner_timers) stacktimer::enable();

    applicator.set_ws_gpu_persistent(ator_ws_gpu_persistent->alloc(applicator.get_wss_gpu_persistent()));
    applicator.preprocess();

    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        data.op_sc->set_ws_persistent(ator_ws_gpu_persistent->alloc(data.op_sc->get_wss_persistent()));
    }

    AllocatorCBMB_new ator_tmp_cbmba(AllocatorGPU_new::get_singleton(), feti.gpu_tmp_mem, feti.gpu_tmp_size);

    ssize_t free_mem_before_preprocess = gpu::mgm::get_device_memory_free();

    stacktimer::push("TotalFETIExplicitGeneralSchurGpu::set preprocess");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        gpu::mgm::queue & q = queues[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        stacktimer::info("TotalFETIExplicitGeneralSchurGpu::set preprocess subdomain %zu", di);

        void * ws_tmp = ator_tmp_cbmba.alloc(data.op_sc->get_wss_tmp_preprocess());

        data.op_sc->preprocess_submit(ws_tmp);

        gpu::mgm::submit_host_function(q, [&,ws_tmp](){
            void * ws_tmp_ = ws_tmp;
            ator_tmp_cbmba.free(ws_tmp_);
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
void TotalFETIExplicitGeneralSchurGpu<T,I>::update(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralSchurGpu::update");

    AllocatorCBMB_new ator_tmp_cbmba(AllocatorGPU_new::get_singleton(), feti.gpu_tmp_mem, feti.gpu_tmp_size);

    stacktimer::push("update_mainloop");
    if(cfg.mainloop_update_split == 'C') {
        stacktimer::push("update_mainloop_combined");
        if(!cfg.inner_timers) stacktimer::disable();
        #pragma omp parallel for schedule(static,1) if(cfg.parallel_update)
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_stuff & data = domain_data[di];
            gpu::mgm::queue & q = queues[di % n_queues];

            math::sumCombined(data.Kreg_old, T{1.0}, feti.K[di], feti.RegMat[di]);

            stacktimer::push("update_cbmba_allocation");
            void * ws_tmp = ator_tmp_cbmba.alloc(data.op_sc->get_wss_tmp_perform());
            stacktimer::pop();

            data.op_sc->perform_1_submit();
            data.op_sc->perform_2_submit(ws_tmp);

            gpu::mgm::submit_host_function(q, [&,ws_tmp](){
                void * ws_tmp_ = ws_tmp;
                ator_tmp_cbmba.free(ws_tmp_);
            });

            applicator.update_F(di);
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
        stacktimer::push("update_mainloop_separate_1");
        if(!cfg.inner_timers) stacktimer::disable();
        #pragma omp parallel for schedule(static,1) if(cfg.parallel_update)
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_stuff & data = domain_data[di];

            math::sumCombined(data.Kreg_old, T{1.0}, feti.K[di], feti.RegMat[di]);

            data.op_sc->perform_1_submit();
        }
        if(!cfg.inner_timers) stacktimer::enable();
        stacktimer::pop();

        // clean up the mess from buggy openmp in clang
        utils::run_dummy_parallel_region();

        stacktimer::push("update_mainloop_separate_2");
        stacktimer::push("update_mainloop_separate_2_submit");
        if(!cfg.inner_timers) stacktimer::disable();
        #pragma omp parallel for schedule(static,1) if(cfg.parallel_update)
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_stuff & data = domain_data[di];
            gpu::mgm::queue & q = queues[di % n_queues];

            stacktimer::push("update_cbmba_allocation");
            void * ws_tmp = ator_tmp_cbmba.alloc(data.op_sc->get_wss_tmp_perform());
            stacktimer::pop();

            domain_data[di].op_sc->perform_2_submit(ws_tmp);

            gpu::mgm::submit_host_function(q, [&,ws_tmp](){
                void * ws_tmp_ = ws_tmp;
                ator_tmp_cbmba.free(ws_tmp_);
            });

            applicator.update_F(di);
        }
        if(!cfg.inner_timers) stacktimer::enable();
        stacktimer::pop();

        if(cfg.gpu_wait_after_mainloop_update) {
            stacktimer::push("update_mainloop_separate_2_wait");
            gpu::mgm::device_wait();
            stacktimer::pop();
        }
    }
    stacktimer::pop();

    // clean up the mess from buggy openmp in clang
    utils::run_dummy_parallel_region();

    stacktimer::push("update_compute_vector_d");
    if(!cfg.inner_timers) stacktimer::disable();
    {
        if (feti.updated.B) {
            d.resize();
        }
        std::vector<Vector_Dense<T,I>> Kplus_fs(n_domains);
        #pragma omp parallel for schedule(static,1) if(cfg.parallel_update)
        for(size_t di = 0; di < n_domains; di++) {
            Kplus_fs[di].resize(feti.f[di].size);
            VectorDenseView_new<T> f_new = VectorDenseView_new<T>::from_old(feti.f[di]);
            VectorDenseView_new<T> Kplus_f_new = VectorDenseView_new<T>::from_old(Kplus_fs[di]);
            domain_data[di].op_sc->solve_A11(f_new, Kplus_f_new);
        }
        applyB(feti, Kplus_fs, d);
        d.synchronize();
        math::add(d, T{-1}, feti.c);
    }
    if(!cfg.inner_timers) stacktimer::enable();
    stacktimer::pop();

    stacktimer::push("update_final_wait");
    gpu::mgm::device_wait();
    stacktimer::pop();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurGpu<T,I>::apply(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralSchurGpu::apply (vector)");

    VectorDenseView_new<T> x_cluster_new = VectorDenseView_new<T>::from_old(x_cluster);
    VectorDenseView_new<T> y_cluster_new = VectorDenseView_new<T>::from_old(y_cluster);

    applicator.apply(x_cluster_new, y_cluster_new, feti.gpu_tmp_mem, feti.gpu_tmp_size);

    y_cluster.synchronize();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurGpu<T,I>::apply(const Matrix_Dual<T> &X_cluster, Matrix_Dual<T> &Y_cluster)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralSchurGpu::apply (matrix)");

    // X_cluster and Y_cluster have the indivudual dual vectors contiguously in rows, but applicator assumes they are contiguously in columns

    MatrixDenseView_new<T> X_cluster_new = MatrixDenseView_new<T>::from_old(X_cluster).get_transposed_reordered_view();
    MatrixDenseView_new<T> Y_cluster_new = MatrixDenseView_new<T>::from_old(Y_cluster).get_transposed_reordered_view();

    applicator.apply(X_cluster_new, Y_cluster_new, feti.gpu_tmp_mem, feti.gpu_tmp_size);
    
    Y_cluster.synchronize();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurGpu<T,I>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralSchurGpu::toPrimal");
    if(!cfg.inner_timers) stacktimer::disable();

    #pragma omp parallel for schedule(static,1)
    for (size_t di = 0; di < n_domains; ++di) {
        Vector_Dense<T,I> z;
        z.resize(y[di]);
        applyBt(feti, di, x, z, T{-1});
        math::add(z, T{1}, feti.f[di]);
        VectorDenseView_new<T> z_new = VectorDenseView_new<T>::from_old(z);
        VectorDenseView_new<T> y_new = VectorDenseView_new<T>::from_old(y[di]);
        domain_data[di].op_sc->solve_A11(z_new, y_new);
    }

    if(!cfg.inner_timers) stacktimer::enable();
    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurGpu<T,I>::print(const step::Step &step)
{
    eslog::error("TotalFETIExplicitGeneralSchurGpu::print not implemented");
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurGpu<T,I>::setup_config(config & cfg, const FETIConfiguration & feti_ecf_config)
{
    // defaults are set in config definition
    // if ecf value is auto, cfg value is not changed

    using ecf_config = DualopTotalfetiExplicitGeneralSchurGpuConfig;
    const ecf_config & ecf = feti_ecf_config.dualop_totalfeti_explicit_generalschur_gpu_config;

    switch(ecf.parallel_set) {
        case ecf_config::AUTOBOOL::AUTO: break;
        case ecf_config::AUTOBOOL::TRUE:  cfg.parallel_set = true;  break;
        case ecf_config::AUTOBOOL::FALSE: cfg.parallel_set = false; break;
    }

    switch(ecf.parallel_update) {
        case ecf_config::AUTOBOOL::AUTO: break;
        case ecf_config::AUTOBOOL::TRUE:  cfg.parallel_update = true;  break;
        case ecf_config::AUTOBOOL::FALSE: cfg.parallel_update = false; break;
    }

    switch(ecf.parallel_apply) {
        case ecf_config::AUTOBOOL::AUTO: break;
        case ecf_config::AUTOBOOL::TRUE:  cfg.parallel_apply = true;  break;
        case ecf_config::AUTOBOOL::FALSE: cfg.parallel_apply = false; break;
    }

    switch(ecf.mainloop_update_split) {
        case ecf_config::MAINLOOP_UPDATE_SPLIT::AUTO: break;
        case ecf_config::MAINLOOP_UPDATE_SPLIT::COMBINED: cfg.mainloop_update_split = 'C'; break;
        case ecf_config::MAINLOOP_UPDATE_SPLIT::SEPARATE: cfg.mainloop_update_split = 'S'; break;
    }

    switch(ecf.synchronize_after_update_mainloop) {
        case ecf_config::AUTOBOOL::AUTO: break;
        case ecf_config::AUTOBOOL::TRUE:  cfg.gpu_wait_after_mainloop_update = true;  break;
        case ecf_config::AUTOBOOL::FALSE: cfg.gpu_wait_after_mainloop_update = false; break;
    }

    switch(ecf.timers_outer) {
        case ecf_config::AUTOBOOL::AUTO: break;
        case ecf_config::AUTOBOOL::TRUE:  cfg.outer_timers = true;  break;
        case ecf_config::AUTOBOOL::FALSE: cfg.outer_timers = false; break;
    }

    switch(ecf.timers_inner) {
        case ecf_config::AUTOBOOL::AUTO: break;
        case ecf_config::AUTOBOOL::TRUE:  cfg.inner_timers = true;  break;
        case ecf_config::AUTOBOOL::FALSE: cfg.inner_timers = false; break;
    }

    switch(ecf.print_config) {
        case ecf_config::AUTOBOOL::AUTO: break;
        case ecf_config::AUTOBOOL::TRUE:  cfg.print_config = true;  break;
        case ecf_config::AUTOBOOL::FALSE: cfg.print_config = false; break;
    }

    switch(ecf.order_F) {
        case ecf_config::MATRIX_ORDER::AUTO: break;
        case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.order_F = 'R'; break;
        case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.order_F = 'C'; break;
    }

    switch(ecf.schur_impl) {
        case ecf_config::SCHUR_IMPL::AUTO: break;
        case ecf_config::SCHUR_IMPL::MANUAL_SIMPLE: cfg.schur_impl = schur_impl_t::manual_simple; break;
        case ecf_config::SCHUR_IMPL::TRIANGULAR:    cfg.schur_impl = schur_impl_t::triangular;    break;
    }

    switch(ecf.apply_where) {
        case ecf_config::CPU_GPU::AUTO: cfg.apply_where = 'G'; break;
        case ecf_config::CPU_GPU::CPU:  cfg.apply_where = 'C'; break;
        case ecf_config::CPU_GPU::GPU:  cfg.apply_where = 'G'; break;
    }
}



#define INSTANTIATE_T_I(T,I) \
template class TotalFETIExplicitGeneralSchurGpu<T,I>;

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

