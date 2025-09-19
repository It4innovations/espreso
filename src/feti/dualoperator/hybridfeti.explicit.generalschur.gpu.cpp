
#include "feti/dualoperator/hybridfeti.explicit.generalschur.gpu.h"

#include "math/primitives_new/allocator_new.h"
#include "feti/common/applyB.h"
#include "basis/utilities/sysutils.h"
#include "basis/utilities/minmaxavg.h"
#include "basis/utilities/stacktimer.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/clusterstore.h"
#include "math/operations/lincomb_vector.h"
#include "math/operations/copy_dnx.h"
#include "math/operations/lincomb_matrix_dnx.h"
#include "math/operations/submatrix_dnx_dnx_noncontig.h"
#include "math/operations/supermatrix_dnx_dnx_noncontig.h"

#include <algorithm>



namespace espreso {



template<typename T, typename I>
HybridFETIExplicitGeneralSchurGpu<T,I>::HybridFETIExplicitGeneralSchurGpu(FETI<T> &feti) : DualOperator<T>(feti)    , main_q(feti.main_q)
    , queues(feti.queues)
    , handles_dense(feti.handles_dense)
    , handles_sparse(feti.handles_sparse)
{
    setup_config(cfg, feti.configuration);
}



template<typename T, typename I>
HybridFETIExplicitGeneralSchurGpu<T,I>::~HybridFETIExplicitGeneralSchurGpu()
{
}



template<typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::info()
{
    eslog::info(" = EXPLICIT HYBRID FETI OPERATOR ON GPU                                                      = \n");

    if(cfg.print_config) {
        auto order_to_string = [](char order){ switch(order){ case 'R': return "ROW_MAJOR"; case 'C': return "COL_MAJOR"; default: return "UNDEFINED"; }};
        auto bool_to_string = [](bool val){ return val ? "TRUE" : "FALSE";};
        auto loop_split_to_string = [](char val){ switch(val){ case 'C': return "COMBINED"; case 'S': return "SEPARATE"; default: return "UNDEFINED"; }};
        auto schur_impl_to_string = [](schur_impl_t schur_impl){ switch(schur_impl) { case schur_impl_t::autoselect: return "autoselect"; case schur_impl_t::manual_simple: return "manual_simple"; case schur_impl_t::triangular: return "triangular"; default: return "UNDEFINED"; } };
        auto schur_impl_to_string_actual = [](schur_impl_t schur_impl){ return gpu::operations::schur_hcsx_ddny<T,I>::make(schur_impl)->get_name(); };

        eslog::info(" =   %-50s       %+30s = \n", "mainloop_update_split", loop_split_to_string(cfg.mainloop_update_split));
        eslog::info(" =   %-50s       %+30s = \n", "synchronize_after_update_mainloops", bool_to_string(cfg.gpu_wait_after_mainloop_update));
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
void HybridFETIExplicitGeneralSchurGpu<T,I>::setup()
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIExplicitGeneralSchurGpu::setup");

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

    d_F1s_allocated.resize((n_domains - 1) / 2 + 1);
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
            MatrixDenseData_new<T> & d_F1_allocd = d_F1s_allocated[allocated_F_index];
            per_domain_stuff & data_bigger = domain_data[di_bigger];
            per_domain_stuff & data_di = domain_data[di];
            if(di == di_bigger) {
                d_F1_allocd.set(data_bigger.n_dofs_interface + (cfg.order_F == 'R'), data_bigger.n_dofs_interface + (cfg.order_F == 'C'), cfg.order_F, ator_ws_gpu_persistent.get());
            }
            data_di.d_F1.set_view(data_di.n_dofs_interface, data_di.n_dofs_interface, d_F1_allocd.ld, d_F1_allocd.order, nullptr, d_F1_allocd.ator);
            data_di.op_sub_F1_from_allocd.set_matrix_src(&d_F1_allocd);
            data_di.op_sub_F1_from_allocd.set_matrix_dst(&data_di.d_F1);
            if(di == di_bigger) {
                data_di.op_sub_F1_from_allocd.set_bounds(0, data_di.d_F1.nrows, 0, data_di.d_F1.ncols);
                data_di.d_F1.prop.uplo = (d_F1_allocd.order == 'R') ? 'U' : 'L';
            }
            else {
                data_di.op_sub_F1_from_allocd.set_bounds((int)(cfg.order_F == 'R'), data_di.d_F1.nrows + (int)(cfg.order_F == 'R'), (int)(cfg.order_F == 'C'), data_di.d_F1.ncols + (int)(cfg.order_F == 'C'));
                data_di.d_F1.prop.uplo = (d_F1_allocd.order == 'R') ? 'L' : 'U';
            }
            data_di.d_F1.prop.symm = MatrixSymmetry_new::hermitian;
        }
    }

    size_t wss_pers_F1s = 0;
    for(const auto & d_F1 : d_F1s_allocated) {
        wss_pers_F1s += d_F1.get_memory_impact();
    }
    total_wss_gpu_persistent += wss_pers_F1s;


    {
        std::vector<MatrixDenseView_new<T>*> Fs_vector(n_domains);
        for(size_t di = 0; di < n_domains; di++) {
            Fs_vector[di] = &domain_data[di].d_F1;
        }

        F1_applicator.set_config(cfg.inner_timers);
        F1_applicator.set_handles(&main_q, &queues, &handles_dense);
        F1_applicator.set_feti(&feti);
        F1_applicator.set_memory('C', 'G');
        F1_applicator.set_D2C_map(&feti.D2C);
        F1_applicator.set_Fs(Fs_vector);
        F1_applicator.set_apply_target(cfg.apply_where);
        F1_applicator.setup();
        total_wss_gpu_persistent += F1_applicator.get_wss_gpu_persistent();
    }

    stacktimer::push("TotalFETIExplicitGeneralSchurGpu::setup op_sc");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1)
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

        data.B1t = MatrixCsxView_new<T,I>::from_old(feti.B1[di]).get_transposed_reordered_view();
        data.Kreg = MatrixCsxView_new<T,I>::from_old(data.Kreg_old);

        data.op_sc = gpu::operations::schur_hcsx_ddny<T,I>::make(cfg.schur_impl);
        data.op_sc->set_handles(q, hs, hd);
        data.op_sc->set_coefficients(-1);
        data.op_sc->set_matrix(&data.Kreg, &data.B1t, nullptr, nullptr);
        data.op_sc->set_sc(&data.d_F1);
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

    stacktimer::info("TotalFETIExplicitGeneralSchurGpu::setup wss_gpu_persistent_F1s %zd", wss_pers_F1s);
    stacktimer::info("TotalFETIExplicitGeneralSchurGpu::setup wss_gpu_persistent_applicator %zd", F1_applicator.get_wss_gpu_persistent());
    stacktimer::info("TotalFETIExplicitGeneralSchurGpu::setup wss_gpu_persistent_op_sc %zu", op_sc_wss_gpu_persistent);
    stacktimer::info("TotalFETIExplicitGeneralSchurGpu::setup wss_gpu_internal_op_sc %zu", op_sc_wss_gpu_internal);
    stacktimer::info("TotalFETIExplicitGeneralSchurGpu::setup max_wss_tmp_preprocess_op_sc %zu", std::max_element(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & l, const per_domain_stuff & r){ return l.op_sc->get_wss_tmp_preprocess() < r.op_sc->get_wss_tmp_preprocess(); })->op_sc->get_wss_tmp_preprocess());
    stacktimer::info("TotalFETIExplicitGeneralSchurGpu::setup max_wss_tmp_perform_op_sc %zu", std::max_element(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & l, const per_domain_stuff & r){ return l.op_sc->get_wss_tmp_perform() < r.op_sc->get_wss_tmp_perform(); })->op_sc->get_wss_tmp_perform());

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::set(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIExplicitGeneralSchurGpu::set");

    Btx.resize(feti.K.size());
    KplusBtx.resize(feti.K.size());
    dKB0.resize(feti.K.size());
    B0mu.resize(feti.K.size());
    hfetiBtx.resize(feti.K.size());

    if(ws_gpu_persistent == nullptr) eslog::error("persistent gpu workspace is not set\n");

    ator_ws_gpu_persistent->set(ws_gpu_persistent, total_wss_gpu_persistent);

    for(auto & d_F1 : d_F1s_allocated) d_F1.alloc();

    if(!cfg.inner_timers) stacktimer::disable();
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & my = domain_data[di];
        my.op_sub_F1_from_allocd.perform();
    }
    if(!cfg.inner_timers) stacktimer::enable();

    F1_applicator.set_ws_gpu_persistent(ator_ws_gpu_persistent->alloc(F1_applicator.get_wss_gpu_persistent()));
    F1_applicator.preprocess();

    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        data.op_sc->set_ws_persistent(ator_ws_gpu_persistent->alloc(data.op_sc->get_wss_persistent()));
    }

    AllocatorCBMB_new ator_tmp_cbmba(AllocatorGPU_new::get_singleton(), feti.gpu_tmp_mem, feti.gpu_tmp_size);

    ssize_t free_mem_before_preprocess = gpu::mgm::get_device_memory_free();

    stacktimer::push("HybridFETIExplicitGeneralSchurGpu::set preprocess");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1)
    for(size_t di = 0; di < n_domains; di++) {
        gpu::mgm::queue & q = queues[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        stacktimer::info("HybridFETIExplicitGeneralSchurGpu::set preprocess subdomain %zu", di);

        void * ws_tmp = ator_tmp_cbmba.alloc(data.op_sc->get_wss_tmp_preprocess());
    
        data.op_sc->preprocess_submit(ws_tmp);

        gpu::mgm::submit_host_function(q, [&,ws_tmp](){
            void * ws_tmp_ = ws_tmp;
            ator_tmp_cbmba.free(ws_tmp_);
        });

        Btx[di].resize(feti.K[di].nrows);
        KplusBtx[di].resize(feti.K[di].nrows);
        B0mu[di].resize(feti.K[di].nrows);
        hfetiBtx[di].resize(feti.K[di].nrows);
        math::set(Btx[di], T{0});
    }
    gpu::mgm::device_wait();
    if(!cfg.inner_timers) stacktimer::enable();
    stacktimer::pop();

    // clean up the mess from buggy openmp in clang
    utils::run_dummy_parallel_region();

    ssize_t free_mem_after_preprocess = gpu::mgm::get_device_memory_free();
    ssize_t allocd_in_preprocess = free_mem_before_preprocess - free_mem_after_preprocess;
    stacktimer::info("TotalFETIExplicitScTriaGpu::set allocd_in_preprocess %zd", allocd_in_preprocess);

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::update(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIExplicitGeneralSchurGpu::update");

    AllocatorCBMB_new ator_tmp_cbmba(AllocatorGPU_new::get_singleton(), feti.gpu_tmp_mem, feti.gpu_tmp_size);

    stacktimer::push("update_mainloop");
    if(cfg.mainloop_update_split == 'C') {
        stacktimer::push("update_mainloop_combined");
        if(!cfg.inner_timers) stacktimer::disable();
        #pragma omp parallel for schedule(static,1)
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

            F1_applicator.update_F(di);
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
        #pragma omp parallel for schedule(static,1)
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
        #pragma omp parallel for schedule(static,1)
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

            F1_applicator.update_F(di);
        }
        if(!cfg.inner_timers) stacktimer::enable();
        stacktimer::pop();

        if(cfg.gpu_wait_after_mainloop_update) {
            stacktimer::push("update_mainloop_separate_2_wait");
            gpu::mgm::device_wait();
            stacktimer::pop();
        }
        stacktimer::pop();
    }
    stacktimer::pop();

    // clean up the mess from buggy openmp in clang
    utils::run_dummy_parallel_region();

    stacktimer::push("hybrid_feti_things");
    {
        int nR1 = 0, nKR1 = 0;
        for (size_t di = 0; di < feti.K.size(); ++di) {
            nR1 += feti.R1[di].nrows;
            nKR1 += feti.KR1[di].nrows;
        }
        if (nR1 == 0 && nKR1 == 0) {
            eslog::error("HYBRID FETI: provide kernels for B0 gluing.\n");
        }
        isRegularK = nR1 == 0;

        _computeB0();
        _computeF0();
        _computeG0();
        _computeS0();
        g.resize(feti.cluster.gl_size);
        beta.resize(G0.nrows);
        mu.resize(feti.cluster.gl_size);

        if (info::ecf->output.print_matrices) {
            eslog::storedata(" STORE: feti/dualop/{B0, F0, G0, S0}\n");
            for (size_t d = 0; d < feti.B1.size(); ++d) {
                math::store(B0[d], utils::filename(utils::debugDirectory(step) + "/feti/dualop", "B0" + std::to_string(d)).c_str());
                math::store(D2C0[d], utils::filename(utils::debugDirectory(step) + "/feti/dualop", "D2C0" + std::to_string(d)).c_str());
            }
            math::store(F0, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "F0").c_str());
            math::store(G0, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "G0").c_str());
            math::store(S0, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "S0").c_str());
        }
    }
    stacktimer::pop();

    stacktimer::push("update_compute_vector_d");
    if(!cfg.inner_timers) stacktimer::disable();
    {
        if (feti.updated.B) {
            d.resize();
        }
        std::vector<Vector_Dense<T,I>> Kplus_fs(n_domains);
        #pragma omp parallel for schedule(static,1)
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
void HybridFETIExplicitGeneralSchurGpu<T,I>::apply(const Vector_Dual<T> &x_cluster_old, Vector_Dual<T> &y_cluster_old)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIExplicitGeneralSchurGpu::apply (vector)");

    VectorDenseView_new<T> x_cluster = VectorDenseView_new<T>::from_old(x_cluster_old);
    VectorDenseView_new<T> y_cluster = VectorDenseView_new<T>::from_old(y_cluster_old);

    Allocator_new * ator_2 = ((cfg.apply_where == 'G') ? (Allocator_new*)AllocatorHostPinned_new::get_singleton() : (Allocator_new*)AllocatorCPU_new::get_singleton());
    VectorDenseData_new<T> x_cluster_2;
    VectorDenseData_new<T> y_cluster_2;
    x_cluster_2.set(x_cluster.size, ator_2);
    y_cluster_2.set(y_cluster.size, ator_2);
    x_cluster_2.alloc();
    y_cluster_2.alloc();

    std::copy_n(x_cluster.vals, x_cluster.size, x_cluster_2.vals);

    F1_applicator.apply(x_cluster_2, y_cluster_2, feti.gpu_tmp_mem, feti.gpu_tmp_size, [&](){_apply_hfeti_stuff(x_cluster_old, y_cluster_old);});

    math::operations::lincomb_vector<T>::do_all(&y_cluster, T{1}, &y_cluster, T{1}, &y_cluster_2);

    y_cluster_old.synchronize();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::apply(const Matrix_Dual<T> &X_cluster_old, Matrix_Dual<T> &Y_cluster_old)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIExplicitGeneralSchurGpu::apply (matrix)");

    MatrixDenseView_new<T> X_cluster = MatrixDenseView_new<T>::from_old(X_cluster_old).get_transposed_reordered_view();
    MatrixDenseView_new<T> Y_cluster = MatrixDenseView_new<T>::from_old(Y_cluster_old).get_transposed_reordered_view();

    Allocator_new * ator_2 = ((cfg.apply_where == 'G') ? (Allocator_new*)AllocatorHostPinned_new::get_singleton() : (Allocator_new*)AllocatorCPU_new::get_singleton());
    MatrixDenseData_new<T> X_cluster_2;
    MatrixDenseData_new<T> Y_cluster_2;
    X_cluster_2.set(X_cluster.nrows, X_cluster.ncols, X_cluster.order, ator_2);
    Y_cluster_2.set(Y_cluster.nrows, Y_cluster.ncols, Y_cluster.order, ator_2);
    X_cluster_2.alloc();
    Y_cluster_2.alloc();

    math::operations::copy_dnx<T>::do_all(&X_cluster, &X_cluster_2);

    F1_applicator.apply(X_cluster_2, Y_cluster_2, feti.gpu_tmp_mem, feti.gpu_tmp_size, [&](){_apply_hfeti_stuff(X_cluster_old, Y_cluster_old);});

    math::operations::lincomb_matrix_dnx<T>::do_all(&Y_cluster, T{1}, &Y_cluster, T{1}, &Y_cluster_2);

    Y_cluster_old.synchronize();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::apply(const Matrix_Dual<T> &X_cluster_old, Matrix_Dual<T> &Y_cluster_old, const std::vector<int> &filter)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIExplicitGeneralSchurGpu::apply (matrix,filter)");

    MatrixDenseView_new<T> X_cluster = MatrixDenseView_new<T>::from_old(X_cluster_old).get_transposed_reordered_view();
    MatrixDenseView_new<T> Y_cluster = MatrixDenseView_new<T>::from_old(Y_cluster_old).get_transposed_reordered_view();

    std::vector<I> filtered_idxs;
    for(size_t i = 0; i < filter.size(); i++) {
        if(filter[i]) {
            filtered_idxs.push_back(i);
        }
    }
    VectorDenseView_new<I> filter_map;
    filter_map.set_view(filtered_idxs.size(), filtered_idxs.data(), AllocatorDummy_new::get_singleton(true,false));

    Allocator_new * ator_2 = ((cfg.apply_where == 'G') ? (Allocator_new*)AllocatorHostPinned_new::get_singleton() : (Allocator_new*)AllocatorCPU_new::get_singleton());
    MatrixDenseData_new<T> X_cluster_2;
    MatrixDenseData_new<T> Y_cluster_2;
    X_cluster_2.set(X_cluster.nrows, X_cluster.ncols, X_cluster.order, ator_2);
    Y_cluster_2.set(Y_cluster.nrows, Y_cluster.ncols, Y_cluster.order, ator_2);
    X_cluster_2.alloc();
    Y_cluster_2.alloc();

    math::operations::submatrix_dnx_dnx_noncontig<T,I>::do_all(&X_cluster, &X_cluster_2, nullptr, &filter_map);

    F1_applicator.apply(X_cluster_2, Y_cluster_2, feti.gpu_tmp_mem, feti.gpu_tmp_size, [&](){_apply_hfeti_stuff(X_cluster_old, Y_cluster_old, filter);});

    math::operations::supermatrix_dnx_dnx_noncontig<T,I>::do_all(&Y_cluster_2, &Y_cluster, nullptr, &filter_map, math::operations::supermatrix_dnx_dnx_noncontig<T,I>::mode::accumulate);

    Y_cluster_old.synchronize();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIExplicitGeneralSchurGpu::toPrimal");
    if(!cfg.inner_timers) stacktimer::disable();

    #pragma omp parallel for schedule(static,1)
    for (size_t di = 0; di < feti.K.size(); ++di) {
        applyBt(feti, di, x, Btx[di]);
        math::copy(KplusBtx[di], feti.f[di]);
        math::add(KplusBtx[di], T{-1}, Btx[di]);
    }
    _applyK(KplusBtx, y, true);

    if(!cfg.inner_timers) stacktimer::enable();
    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template <typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::BtL(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIExplicitGeneralSchurGpu::BtL");

    #pragma omp parallel for schedule(static,1)
    for (size_t di = 0; di < feti.K.size(); ++di) {
        applyBt(feti, di, x, y[di]);
    }

    _compute_beta_mu(y);

    #pragma omp parallel for schedule(static,1)
    for (size_t di = 0; di < feti.K.size(); ++di) {
        // y = (b - B0t * mu)
        math::spblas::applyT(y[di], T{-1}, B0[di], D2C0[di].data(), mu);
    }

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template <typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::_apply_hfeti_stuff(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    for (size_t di = 0; di < feti.K.size(); ++di) {
        applyBt(feti, di, x, Btx[di]);
    }

    _applyK(Btx, KplusBtx, false);
    applyB(feti, KplusBtx, y);
}



template <typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::_apply_hfeti_stuff(const Matrix_Dual<T> &x, Matrix_Dual<T> &y)
{
    Vector_Dual<T> _x, _y;
    _x.size = _y.size = x.ncols;

    for (int r = 0; r < x.nrows; ++r) {
        _x.vals = x.vals + x.ncols * r;
        _y.vals = y.vals + y.ncols * r;

        for (size_t di = 0; di < feti.K.size(); ++di) {
            applyBt(feti, di, _x, Btx[di]);
        }

        _applyK(Btx, KplusBtx, false);
        applyB(feti, KplusBtx, _y);
    }
}



template <typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::_apply_hfeti_stuff(const Matrix_Dual<T> &x, Matrix_Dual<T> &y, const std::vector<int> &filter)
{
    Vector_Dual<T> _x, _y;
    _x.size = _y.size = x.ncols;

    for (int r = 0; r < x.nrows; ++r) {
        if (filter[r]) {
            _x.vals = x.vals + x.ncols * r;
            _y.vals = y.vals + y.ncols * r;

            for (size_t di = 0; di < feti.K.size(); ++di) {
                applyBt(feti, di, _x, Btx[di]);
            }

            _applyK(Btx, KplusBtx, false);
            applyB(feti, KplusBtx, _y);
        }
    }
}



// https://dl.acm.org/doi/pdf/10.1145/2929908.2929909

template <typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::_applyK(std::vector<Vector_Dense<T> > &b, std::vector<Vector_Dense<T> > &x, bool do_Kplus_solve)
{
    if(do_Kplus_solve) {
        // x = (K+)^-1 * b
        // actually x = K+ * b = (Kreg)^-1 * b
        #pragma omp parallel for schedule(static,1)
        for (size_t di = 0; di < feti.K.size(); ++di) {
            // KSolver[di].solve(b[di], x[di]);
            x[di].resize(b[di].size);
            VectorDenseView_new<T> rhs = VectorDenseView_new<T>::from_old(b[di]);
            VectorDenseView_new<T> sol = VectorDenseView_new<T>::from_old(x[di]);
            domain_data[di].op_sc->solve_A11(rhs, sol);
        }
    }
    else {
        #pragma omp parallel for  schedule(static,1)
        for (size_t di = 0; di < feti.K.size(); ++di) {
            std::fill_n(x[di].vals, x[di].size, T{0});
        }
    }

    _compute_beta_mu(b);

    // x -= (K+)^-1 * (B0t * mu)
    #pragma omp parallel for schedule(static,1)
    for (size_t di = 0; di < feti.K.size(); ++di) {
        Vector_Dense<T> mu_di; mu_di.resize(D2C0[di].size());
        for (size_t i = 0; i < D2C0[di].size(); ++i) {
            mu_di.vals[i] = mu.vals[D2C0[di][i]];
        }

        math::blas::multiply(T{1}, dKB0[di], mu_di, T{0}, hfetiBtx[di], true);
        math::add(x[di], T{-1}, hfetiBtx[di]);
    }

    // x += R * beta
    if (!isRegularK) {
        #pragma omp parallel for schedule(static,1)
        for (size_t di = 0; di < feti.K.size(); ++di) {
            Vector_Dense<T> _beta;
            _beta.size = origR1[di].nrows;
            _beta.vals = beta.vals + G0offset[di];
            math::blas::multiply(T{1}, origR1[di], _beta, T{1}, x[di], true);
        }
    }
}



template <typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::_compute_beta_mu(std::vector<Vector_Dense<T> > &b)
{
    // g = B0 * (K+)^-1 * b
    math::set(g, T{0});
    for (size_t di = 0; di < feti.K.size(); ++di) {
        Vector_Dense<T> KB0b; KB0b.resize(dKB0[di].nrows);
        math::blas::apply(KB0b, T{1}, dKB0[di], T{0}, b[di]);
        for (size_t i = 0; i < D2C0[di].size(); ++i) {
            g.vals[D2C0[di][i]] += KB0b.vals[i];
        }

        if (!isRegularK) {
            // e = Rt * b
            Vector_Dense<T> _e;
            _e.size = origR1[di].nrows;
            _e.vals = beta.vals + G0offset[di];
            math::blas::multiply(T{1}, origR1[di], b[di], T{0}, _e); // assemble -e
        }
    }

    if (!isRegularK) {
        auto &F0g = mu; // mu is tmp variable that can be used here
        F0Solver.solve(g, F0g);
        // e = G0 * F^-1 * g - e
        math::spblas::apply(beta, T{1}, G0, F0g);
        // beta = (S+)^-1 * (G0 * F0^-1 * g - e)
        Splus.solve(beta);
        // g = g - G0t * beta
        math::spblas::applyT(g, T{-1}, G0, beta);
    }
    // mu = F0^-1 * (g - G0t * beta)
    F0Solver.solve(g, mu);
}



template <typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::_computeB0()
{
if constexpr(utils::is_real<T>()) {
    std::vector<Matrix_Dense<T> > &R = isRegularK ? feti.KR1 : feti.R1;

    B0.clear();
    B0.resize(feti.K.size());
    D2C0.clear();
    D2C0.resize(feti.K.size());

    struct __ijv__ {
        int i,j; double v;

        bool operator<(const __ijv__ &other) const {
            if (i == other.i) {
                return j < other.j;
            }
            return i < other.i;
        }
    };

    std::vector<int> rindex;
    std::vector<std::vector<__ijv__> > iB0(B0.size());

    auto dual = info::mesh->domains->localDual->begin();
    std::vector<int> rows(info::mesh->clusters->size);
    for (int d1 = 0; d1 < info::mesh->domains->size; ++d1, ++dual) {
        int cluster = info::mesh->domains->cluster[d1];
        for (auto dit = dual->begin(); dit != dual->end(); ++dit) {
            if (d1 < *dit) {
                rindex.push_back(rows[cluster]);
                if (R[d1].nrows < R[*dit].nrows) {
                    rows[cluster] += R[*dit].nrows;
                } else {
                    rows[cluster] += R[d1].nrows ? R[d1].nrows : 1;
                }
            } else {
                auto dualbegin = info::mesh->domains->localDual->begin();
                auto it = std::lower_bound((dualbegin + *dit)->begin(), (dualbegin + *dit)->end(), d1);
                rindex.push_back(rindex[it - dualbegin->begin()]);
            }
        }
    }
    feti.cluster.gl_size = rows.front(); // TODO: more clusters

    dual = info::mesh->domains->localDual->begin();
    for (auto dmap = feti.decomposition->dmap->cbegin(); dmap != feti.decomposition->dmap->cend(); ++dmap) {
        for (auto di1 = dmap->begin(); di1 != dmap->end(); ++di1) {
            for (auto di2 = di1 + 1; di2 != dmap->end(); ++di2) {
                if (feti.decomposition->ismy(di1->domain) && feti.decomposition->ismy(di2->domain)) {
                    auto it = std::lower_bound(
                            (dual + (di1->domain - feti.decomposition->dbegin))->begin(),
                            (dual + (di1->domain - feti.decomposition->dbegin))->end(),
                            di2->domain - feti.decomposition->dbegin);

                    if (it != (dual + (di1->domain - feti.decomposition->dbegin))->end() && *it == di2->domain - feti.decomposition->dbegin) {
                        int d1, d2, d1index, d2index;
                        if (di1->domain < di2->domain) {
                            d1 = di1->domain - feti.decomposition->dbegin;
                            d2 = di2->domain - feti.decomposition->dbegin;
                            d1index = di1->index;
                            d2index = di2->index;
                        } else {
                            d1 = di2->domain - feti.decomposition->dbegin;
                            d2 = di1->domain - feti.decomposition->dbegin;
                            d1index = di2->index;
                            d2index = di1->index;
                        }
                        if (R[d1].nrows) {
                            for (int r = 0; r < R[d1].nrows; ++r) {
                                iB0[d1].push_back({rindex[it - dual->begin()] + r, d1index,  R[d1].vals[R[d1].ncols * r + d1index]});
                                iB0[d2].push_back({rindex[it - dual->begin()] + r, d2index, -R[d1].vals[R[d1].ncols * r + d1index]});
                            }
                        } else {
                            iB0[d1].push_back({rindex[it - dual->begin()], d1index,  (double)(d1 + 1) / info::mesh->domains->size});
                            iB0[d1].push_back({rindex[it - dual->begin()], d2index, -(double)(d1 + 1) / info::mesh->domains->size});
                        }
                    }
                }
            }
        }
    }

    #pragma omp parallel for schedule(static,1)
    for (int d = 0; d < info::mesh->domains->size; ++d) {
        std::sort(iB0[d].begin(), iB0[d].end());
        if (iB0[d].size()) {
            D2C0[d].push_back(iB0[d][0].i);
        }
        for (size_t i = 1; i < iB0[d].size(); ++i) {
            if (D2C0[d].back() != iB0[d][i].i) {
                D2C0[d].push_back(iB0[d][i].i);
            }
        }
        B0[d].resize(D2C0[d].size(), feti.K[d].ncols, iB0[d].size());
        B0[d].rows[0] = 0;
        B0[d].cols[0] = iB0[d][0].j;
        B0[d].vals[0] = iB0[d][0].v;
        for (size_t i = 1, r = 1; i < iB0[d].size(); ++i) {
            B0[d].cols[i] = iB0[d][i].j;
            B0[d].vals[i] = iB0[d][i].v;
            if (iB0[d][i - 1].i != iB0[d][i].i) {
                B0[d].rows[r++] = i;
            }
        }
        B0[d].rows[D2C0[d].size()] = iB0[d].size();
    }
} else {
    eslog::error("complex not supported now, todo\n");
}
}



template <typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::_computeF0()
{
if constexpr(utils::is_real<T>()) {
    if (F0.type == Matrix_Type::UNSET_INVALID_NONE) {
        std::vector<std::vector<int> > csr(feti.cluster.gl_size);
        for (size_t di = 0; di < feti.K.size(); ++di) {
            for (size_t i = 0; i < D2C0[di].size(); ++i) {
                for (size_t j = i; j < D2C0[di].size(); ++j) {
                    csr[D2C0[di][i]].push_back(D2C0[di][j]);
                }
            }
            dKB0[di].resize(B0[di].nrows, B0[di].ncols);
        }

        #pragma omp parallel for schedule(static,1)
        for (size_t i = 0; i < csr.size(); ++i) {
            utils::sortAndRemoveDuplicates(csr[i]);
        }
        int F0nnz = 0;
        for (size_t i = 0; i < csr.size(); ++i) {
            F0nnz += csr[i].size();
        }
        F0.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
        F0.shape = Matrix_Shape::UPPER;
        F0.resize(feti.cluster.gl_size, feti.cluster.gl_size, F0nnz);
        F0.rows[0] = Indexing::CSR;
        for (size_t i = 0; i < csr.size(); ++i) {
            F0.rows[i + 1] = F0.rows[i] + std::max(csr[i].size(), 1UL);
            for (size_t j = 0; j < csr[i].size(); ++j) {
                F0.cols[F0.rows[i] - Indexing::CSR + j] = csr[i][j] + Indexing::CSR;
            }
        }

        permutationF0.resize(feti.K.size());
        for (size_t di = 0, ri = feti.cluster.gl_size; di < feti.K.size(); ++di, ++ri) {
            for (size_t i = 0; i < D2C0[di].size(); ++i) {
                for (size_t j = i; j < D2C0[di].size(); ++j) {
                    int c = std::lower_bound(csr[D2C0[di][i]].begin(), csr[D2C0[di][i]].end(), D2C0[di][j]) - csr[D2C0[di][i]].begin();
                    permutationF0[di].push_back(F0.rows[D2C0[di][i]] + c - Indexing::CSR);
                }
            }
        }

        eslog::checkpointln("FETI: SET F0");
        F0Solver.commit(F0);
        F0Solver.symbolicFactorization();
        eslog::checkpointln("FETI: F0 SYMBOLIC FACTORIZATION");
    }

    math::set(F0, T{0});
    std::vector<int> pi(feti.K.size());
    #pragma omp parallel for schedule(static,1)
    for (size_t di = 0; di < feti.K.size(); ++di) {
        Matrix_Dense<T> dB0; dB0.resize(B0[di].nrows, B0[di].ncols);
        Matrix_Dense<T> dF0; dF0.resize(B0[di].nrows, B0[di].nrows);
        math::set(dB0, T{0});
        for (size_t i = 0; i < D2C0[di].size(); ++i) {
            for (int c = B0[di].rows[i]; c < B0[di].rows[i + 1]; ++c) {
                dB0.vals[i * dB0.ncols + B0[di].cols[c]] = B0[di].vals[c];
            }
        }
        // dF0 = B0 * (K+)^-1 * B0t // dense due to B0 is local
        // KSolver[di].solve(dB0, dKB0[di]);
        {
            dKB0[di].resize(dB0.nrows, dB0.ncols);
            MatrixDenseView_new<T> rhs = MatrixDenseView_new<T>::from_old(dB0).get_transposed_reordered_view();
            MatrixDenseView_new<T> sol = MatrixDenseView_new<T>::from_old(dKB0[di]).get_transposed_reordered_view();
            domain_data[di].op_sc->solve_A11(rhs, sol);
        }
        math::blas::multiply(T{1}, dB0, dKB0[di], T{0}, dF0, false, true);
        // dF0 to F0
        for (size_t i = 0, k = 0; i < D2C0[di].size(); k += ++i) {
            for (size_t j = i; j < D2C0[di].size(); ++j, ++k) {
                #pragma omp atomic
                F0.vals[permutationF0[di][pi[di]++]] += dF0.vals[k];
            }
        }
    }
    eslog::checkpointln("FETI: UPDATE F0");
    F0Solver.numericalFactorization();
    eslog::checkpointln("FETI: F0 NUMERIC FACTORIZATION");
} else {
    eslog::error("complex not supported now, todo\n");
}
}



template <typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::_computeG0()
{
    if (isRegularK) {
        G0.resize(0, feti.cluster.gl_size, 0);
        G0.rows[0] = 0;
        return;
    }

    int G0nnz = 0, G0rows = 0;
    for (size_t di = 0; di < feti.K.size(); ++di) {
        G0nnz += feti.R1[di].nrows * D2C0[di].size();
        G0rows += feti.R1[di].nrows;
    }
    G0.resize(G0rows, feti.cluster.gl_size, G0nnz);
    G0.rows[0] = 0;
    for (size_t di = 0, ri = 0, ci = 0; di < feti.K.size(); ++di) {
        G0offset.push_back(ri);
        for (int r = 0; r < feti.R1[di].nrows; ++r, ++ri) {
            for (size_t c = 0; c < D2C0[di].size(); ++c) {
                G0.cols[ci++] = D2C0[di][c];
            }
            G0.rows[ri + 1] = ci;
        }
    }

    #pragma omp parallel for schedule(static,1)
    for (size_t di = 0; di < feti.K.size(); ++di) {
        // G0 = R1 * B0t
        // dense version: math::blas::multiply(T{1}, feti.R1[di], dB0[t], T{0}, dKB0[t], false, true);
        // sparse version
        for (int r = 0, ri = G0offset[di]; r < feti.R1[di].nrows; ++r, ++ri) {
            for (size_t c = 0; c < D2C0[di].size(); ++c) {
                G0.vals[G0.rows[ri] + c] = 0;
                for (int k = B0[di].rows[c]; k < B0[di].rows[c + 1]; ++k) {
                    G0.vals[G0.rows[ri] + c] -= feti.R1[di].vals[r * feti.R1[di].ncols + B0[di].cols[k]] * B0[di].vals[k];
                }
            }
        }
    }
    eslog::checkpointln("FETI: UPDATE G0");
}



template <typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::_computeS0()
{
    if (isRegularK) {
        S0.resize(0, 0);
        return;
    }

    Matrix_Dense<T> dG0, dF0G0;
    math::copy(dG0, G0);
    F0Solver.solve(dG0, dF0G0);
    S0.resize(G0.nrows, G0.nrows); S0.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    math::blas::multiply(T{1}, dG0, dF0G0, T{0}, S0, false, true); // TODO: use sparse G0?

    origR1.resize(feti.K.size());
    // make S0 regular
    switch (feti.configuration.regularization) {
    case FETIConfiguration::REGULARIZATION::ANALYTIC: {
        T avg = 0;
        for (int r = 0; r < S0.nrows; ++r) {
            avg += S0.vals[S0.nrows * r + r];
        }
        avg /= S0.nrows;

        // kernels: [3, 2, 3]
        // 1 0 0   *  1 0 0 1 0 1 0 0  =   1 0 0 1 0 1 0 0
        // 0 1 0      0 1 0 0 1 0 1 0      0 1 0 0 1 0 1 0
        // 0 0 1      0 0 1 0 0 0 0 1      0 0 1 0 0 0 0 1
        // 1 0 0                           1 0 0 1 0 1 0 0
        // 0 1 0                           0 1 0 0 1 0 1 0
        // 1 0 0                           1 0 0 1 0 1 0 0
        // 0 1 0                           0 1 0 0 1 0 1 0
        // 0 0 1                           0 0 1 0 0 0 0 1

        for (size_t di = 0, ri = 0; di < feti.K.size(); ++di) {
            origR1[di].shallowCopy(feti.R1[di]);
            for (int r = 0; r < feti.R1[di].nrows; ++r, ++ri) {
                for (size_t dj = 0; dj < feti.K.size(); ++dj) {
                    if (r < feti.R1[dj].nrows) {
                        S0.vals[ri * S0.ncols + G0offset[dj] + r] += 0.5 * avg;
                    }
                }
            }
        }
    } break;
    case FETIConfiguration::REGULARIZATION::ALGEBRAIC:
    case FETIConfiguration::REGULARIZATION::SVD: {
        int maxDefect = 0;
        for (size_t di = 0; di < feti.K.size(); ++di) {
            maxDefect = std::max(maxDefect, feti.R1[di].nrows);
            if (maxDefect != feti.R1[di].nrows) {
                eslog::error("Some domains are regular\n");
            }
        }
        Matrix_Dense<T> R;
        Matrix_IJV<T> regMat;
        math::getKernel<T, int>(S0, R, regMat, maxDefect, feti.configuration.sc_size);
        for (int i = 0; i < regMat.nnz; ++i) {
            S0.vals[(regMat.rows[i] - Indexing::IJV) * S0.ncols + regMat.cols[i] - Indexing::IJV] += regMat.vals[i];
        }
        for (size_t di = 0; di < feti.K.size(); ++di) {
            if (feti.configuration.projector_opt & FETIConfiguration::PROJECTOR_OPT::FULL) {
                origR1[di].shallowCopy(feti.R1[di]);
            } else {
                // orthogonalize R for HFETI projector
                origR1[di].resize(feti.R1[di].nrows, feti.R1[di].ncols);
                math::copy(origR1[di], feti.R1[di]);
                Matrix_Dense<T> _R;
                math::lapack::submatrix<T, int>(R, _R, 0, R.nrows, di * R.nrows, (di + 1) * R.nrows);
                math::blas::multiply(T{1}, _R, origR1[di], T{0}, feti.R1[di]);
            }
        }
    } break;
    }

    eslog::checkpointln("FETI: UPDATE S0");

    Splus.commit(S0);
    Splus.factorization();
    eslog::checkpointln("FETI: S0 FACTORIZATION");
}



template<typename T, typename I>
void HybridFETIExplicitGeneralSchurGpu<T,I>::setup_config(config & cfg, const FETIConfiguration & feti_ecf_config)
{
    // defaults are set in config definition
    // if ecf value is auto, cfg value is not changed

    using ecf_config = DualopHybridfetiExplicitGeneralSchurGpuConfig;
    const ecf_config & ecf = feti_ecf_config.dualop_hybridfeti_explicit_generalschur_gpu_config;

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
template class HybridFETIExplicitGeneralSchurGpu<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        /* INSTANTIATE_T(std::complex<double>) */

    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I

}

