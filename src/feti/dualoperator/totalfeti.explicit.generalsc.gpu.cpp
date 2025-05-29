
#include "feti/dualoperator/totalfeti.explicit.generalsc.gpu.h"

#include "feti/common/applyB.h"
#include "basis/utilities/minmaxavg.h"
#include "basis/utilities/stacktimer.h"
#include "esinfo/meshinfo.h"
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

static void replace_if_default(char & param, char deflt)
{
    if(param == '_') {
        param = deflt;
    }
}

template<typename T, typename I>
static void replace_unset_configs(typename TotalFETIExplicitGeneralScGpu<T,I>::config & cfg_dualop)
{
    replace_if_default(cfg_dualop.order_F, 'R');
    replace_if_default(cfg_dualop.mainloop_update_split, 'C');
}

template<typename T, typename I>
static void setup_configs(typename TotalFETIExplicitGeneralScGpu<T,I>::config & cfg_dualop)
{
    set_by_env(cfg_dualop.order_F,                        "ESPRESO_CONFIG_DUALOP_GENERALSC_order_F");
    set_by_env(cfg_dualop.parallel_set,                   "ESPRESO_CONFIG_DUALOP_GENERALSC_parallel_set");
    set_by_env(cfg_dualop.parallel_update,                "ESPRESO_CONFIG_DUALOP_GENERALSC_parallel_update");
    set_by_env(cfg_dualop.parallel_apply,                 "ESPRESO_CONFIG_DUALOP_GENERALSC_parallel_apply");
    set_by_env(cfg_dualop.mainloop_update_split,          "ESPRESO_CONFIG_DUALOP_GENERALSC_mainloop_update_split");
    set_by_env(cfg_dualop.gpu_wait_after_mainloop_update, "ESPRESO_CONFIG_DUALOP_GENERALSC_gpu_wait_after_mainloop_update");
    set_by_env(cfg_dualop.inner_timers,                   "ESPRESO_CONFIG_DUALOP_GENERALSC_inner_timers");
    set_by_env(cfg_dualop.outer_timers,                   "ESPRESO_CONFIG_DUALOP_GENERALSC_outer_timers");
    set_by_env(cfg_dualop.print_parameters,               "ESPRESO_CONFIG_DUALOP_GENERALSC_print_parameters");
    cfg_dualop.sc_is = TotalFETIExplicitGeneralScGpu<T,I>::sc_is_t::autoselect;
}



template<typename T, typename I>
TotalFETIExplicitGeneralScGpu<T,I>::TotalFETIExplicitGeneralScGpu(FETI<T> &feti) : DualOperator<T>(feti)
    , main_q(feti.main_q)
    , queues(feti.queues)
    , handles_dense(feti.handles_dense)
    , handles_sparse(feti.handles_sparse)
{
    setup_configs<T,I>(cfg);
}



template<typename T, typename I>
TotalFETIExplicitGeneralScGpu<T,I>::~TotalFETIExplicitGeneralScGpu()
{
}



template<typename T, typename I>
void TotalFETIExplicitGeneralScGpu<T,I>::info()
{
    eslog::info(" = EXPLICIT TOTAL FETI OPERATOR USING GENERAL SC ON GPU                               = \n");

    if(cfg.print_parameters) {
        auto order_to_string = [](char order){ switch(order){ case 'R': return "ROW_MAJOR"; case 'C': return "COL_MAJOR"; default: return "UNDEFINED"; }};
        auto bool_to_string = [](bool val){ return val ? "TRUE" : "FALSE";};
        auto loop_split_to_string = [](char val){ switch(val){ case 'C': return "COMBINED"; case 'S': return "SEPARATE"; default: return "UNDEFINED"; }};

        eslog::info(" =   %-50s       %+30s = \n", "parallel_set", bool_to_string(cfg.parallel_set));
        eslog::info(" =   %-50s       %+30s = \n", "parallel_set", bool_to_string(cfg.parallel_update));
        eslog::info(" =   %-50s       %+30s = \n", "parallel_set", bool_to_string(cfg.parallel_apply));
        eslog::info(" =   %-50s       %+30s = \n", "mainloop_update_split", loop_split_to_string(cfg.mainloop_update_split));
        eslog::info(" =   %-50s       %+30s = \n", "order_F", order_to_string(cfg.order_F));
    }

    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_domain; }).to_string("  Domain volume [dofs]").c_str());
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_interface; }).to_string("  Domain surface [dofs]").c_str());
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}



template<typename T, typename I>
void TotalFETIExplicitGeneralScGpu<T,I>::set(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralScGpu::set");

    n_domains = feti.K.size();
    n_queues = queues.size();

    domain_data.resize(n_domains);

    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        data.n_dofs_domain = feti.B1[di].ncols;
        data.n_dofs_interface = feti.B1[di].nrows;
    }

    replace_unset_configs<T,I>(cfg);

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
                data_di.d_F = d_F_allocd.get_submatrix_view(0, data_di.n_dofs_interface, 0, data_di.n_dofs_interface);
                data_di.d_F.prop.uplo = 'U';
            }
            data_di.d_F.prop.symm = MatrixSymmetry_new::hermitian;
            if(data_di.d_F.order == 'R') {
                data_di.d_F_old = MatrixDenseView_new<T>::template to_old<I,gpu::mgm::Ad>(data_di.d_F);
            }
            if(data_di.d_F.order == 'C') {
                MatrixDenseView_new<T> F_reordered = data_di.d_F.get_transposed_reordered_view();
                data_di.d_F_old = MatrixDenseView_new<T>::template to_old<I,gpu::mgm::Ad>(F_reordered);
            }
        }
    }
    ssize_t free_mem_after_Falloc = gpu::mgm::get_device_memory_free();
    ssize_t allocd_in_Falloc = free_mem_before_Falloc - free_mem_after_Falloc;
    stacktimer::info("TotalFETIExplicitGeneralScGpu::set allocd_in_Falloc %zd", allocd_in_Falloc);

    {
        std::vector<MatrixDenseView_new<T>*> Fs_vector(n_domains);
        for(size_t di = 0; di < n_domains; di++) {
            Fs_vector[di] = &domain_data[di].d_F;
        }

        applicator.set_handles(&main_q, &queues, &handles_dense);
        applicator.set_dimensions(feti);
        applicator.set_vector_memory('C');
        applicator.set_D2C_map(&feti.D2C);
        applicator.set_Fs(Fs_vector);
        applicator.set_apply_target('G');
        applicator.setup();
        applicator.preprocess();
    }

    stacktimer::push("TotalFETIExplicitGeneralScGpu::set setup");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        gpu::mgm::queue & q = queues[di % n_queues];
        gpu::spblas::handle & hs = handles_sparse[di % n_queues];
        gpu::dnblas::handle & hd = handles_dense[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        stacktimer::info("TotalFETIExplicitGeneralScGpu::set setup subdomain %zu", di);

        math::combine(data.Kreg_old, feti.K[di], feti.RegMat[di]);
        if constexpr(utils::is_real<T>())    data.Kreg_old.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
        if constexpr(utils::is_complex<T>()) data.Kreg_old.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE;
        data.Kreg_old.shape = feti.K[di].shape;

        data.Bt = MatrixCsxView_new<T,I>::from_old(feti.B1[di]).get_transposed_reordered_view();
        data.Kreg = MatrixCsxView_new<T,I>::from_old(data.Kreg_old);

        data.op_sc = gpu::operations::sc_hcsx_ddny<T,I>::make(cfg.sc_is);
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

    total_wss_internal = 0;
    total_wss_persistent = 0;
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        total_wss_internal += data.op_sc->get_wss_internal();
        total_wss_persistent += data.op_sc->get_wss_persistent();
    }
    total_wss_internal = utils::round_up((total_wss_internal * 105) / 100, gpu::mgm::get_natural_pitch_align());

    stacktimer::info("TotalFETIExplicitGeneralScGpu::set total_wss_internal %zu", total_wss_internal);
    stacktimer::info("TotalFETIExplicitGeneralScGpu::set total_wss_persistent %zu", total_wss_persistent);
    stacktimer::info("TotalFETIExplicitGeneralScGpu::set max_wss_tmp_preprocess %zu", std::max_element(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & l, const per_domain_stuff & r){ return l.op_sc->get_wss_tmp_preprocess() < r.op_sc->get_wss_tmp_preprocess(); })->op_sc->get_wss_tmp_preprocess());
    stacktimer::info("TotalFETIExplicitGeneralScGpu::set max_wss_tmp_perform %zu", std::max_element(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & l, const per_domain_stuff & r){ return l.op_sc->get_wss_tmp_perform() < r.op_sc->get_wss_tmp_perform(); })->op_sc->get_wss_tmp_perform());

    ws_persistent = gpu::mgm::memalloc_device(total_wss_persistent);
    ator_ws_persistent = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());
    ator_ws_persistent->set(ws_persistent, total_wss_persistent);

    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        data.op_sc->set_ws_persistent(ator_ws_persistent->alloc(data.op_sc->get_wss_persistent()));
    }

    size_t free_mem = gpu::mgm::get_device_memory_free();
    size_t mem_capacity = gpu::mgm::get_device_memory_capacity();
    size_t reserve = (mem_capacity * 5) / 100;
    wss_tmp_for_cbmba = utils::round_down(free_mem - reserve - total_wss_internal, gpu::mgm::get_natural_pitch_align());
    ws_tmp_for_cbmba = gpu::mgm::memalloc_device(wss_tmp_for_cbmba);
    ator_tmp_cbmba = std::make_unique<AllocatorCBMB_new>(AllocatorGPU_new::get_singleton(), ws_tmp_for_cbmba, wss_tmp_for_cbmba);

    stacktimer::info("TotalFETIExplicitGeneralScGpu::set cbmba_capacity %zu", wss_tmp_for_cbmba);

    ssize_t free_mem_before_preprocess = gpu::mgm::get_device_memory_free();

    stacktimer::push("TotalFETIExplicitGeneralScGpu::set preprocess");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        gpu::mgm::queue & q = queues[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        stacktimer::info("TotalFETIExplicitGeneralScGpu::set preprocess subdomain %zu", di);

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
void TotalFETIExplicitGeneralScGpu<T,I>::update(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralScGpu::update");

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
            void * ws_tmp = ator_tmp_cbmba->alloc(data.op_sc->get_wss_tmp_perform());
            stacktimer::pop();

            data.op_sc->perform_1_submit();
            data.op_sc->perform_2_submit(ws_tmp);

            gpu::mgm::submit_host_function(q, [&,ws_tmp](){
                void * ws_tmp_ = ws_tmp;
                ator_tmp_cbmba->free(ws_tmp_);
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
            void * ws_tmp = ator_tmp_cbmba->alloc(data.op_sc->get_wss_tmp_perform());
            stacktimer::pop();

            domain_data[di].op_sc->perform_2_submit(ws_tmp);

            gpu::mgm::submit_host_function(q, [&,ws_tmp](){
                void * ws_tmp_ = ws_tmp;
                ator_tmp_cbmba->free(ws_tmp_);
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
void TotalFETIExplicitGeneralScGpu<T,I>::_apply(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    VectorDenseView_new<T> x_cluster_new = VectorDenseView_new<T>::from_old(x_cluster);
    VectorDenseView_new<T> y_cluster_new = VectorDenseView_new<T>::from_old(y_cluster);

    applicator.apply(x_cluster_new, y_cluster_new);
}



template<typename T, typename I>
void TotalFETIExplicitGeneralScGpu<T,I>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    _apply(x, y);
    y.synchronize();
}



template<typename T, typename I>
void TotalFETIExplicitGeneralScGpu<T,I>::apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y)
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
void TotalFETIExplicitGeneralScGpu<T,I>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralScGpu::toPrimal");
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
void TotalFETIExplicitGeneralScGpu<T,I>::print(const step::Step &step)
{
    eslog::error("TotalFETIExplicitGeneralScGpu::print not implemented");
}



#define INSTANTIATE_T_I(T,I) \
template class TotalFETIExplicitGeneralScGpu<T,I>;

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

