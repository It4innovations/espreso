
#include "feti/dualoperator/totalfeti.explicit.generalsc.cpu.h"

#include "feti/common/applyB.h"
#include "basis/utilities/minmaxavg.h"
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

template<typename T, typename I>
static void replace_unset_configs(typename TotalFETIExplicitGeneralScCpu<T,I>::config & cfg_dualop)
{
    replace_if_default(cfg_dualop.order_F, 'R');
    replace_if_default(cfg_dualop.mainloop_update_split, 'C');
}

template<typename T, typename I>
static void setup_configs(typename TotalFETIExplicitGeneralScCpu<T,I>::config & cfg_dualop)
{
    set_by_env(cfg_dualop.order_F,                        "ESPRESO_CONFIG_DUALOP_GENERALSC_order_F");
    set_by_env(cfg_dualop.parallel_set,                   "ESPRESO_CONFIG_DUALOP_GENERALSC_parallel_set");
    set_by_env(cfg_dualop.parallel_update,                "ESPRESO_CONFIG_DUALOP_GENERALSC_parallel_update");
    set_by_env(cfg_dualop.parallel_apply,                 "ESPRESO_CONFIG_DUALOP_GENERALSC_parallel_apply");
    set_by_env(cfg_dualop.mainloop_update_split,          "ESPRESO_CONFIG_DUALOP_GENERALSC_mainloop_update_split");
    // set_by_env(cfg_dualop.gpu_wait_after_mainloop_update, "ESPRESO_CONFIG_DUALOP_GENERALSC_gpu_wait_after_mainloop_update");
    set_by_env(cfg_dualop.inner_timers,                   "ESPRESO_CONFIG_DUALOP_GENERALSC_inner_timers");
    set_by_env(cfg_dualop.outer_timers,                   "ESPRESO_CONFIG_DUALOP_GENERALSC_outer_timers");
    set_by_env(cfg_dualop.print_parameters,               "ESPRESO_CONFIG_DUALOP_GENERALSC_print_parameters");
}



template<typename T, typename I>
TotalFETIExplicitGeneralScCpu<T,I>::TotalFETIExplicitGeneralScCpu(FETI<T> &feti) : DualOperator<T>(feti)
{
    setup_configs<T,I>(cfg);
}



template<typename T, typename I>
TotalFETIExplicitGeneralScCpu<T,I>::~TotalFETIExplicitGeneralScCpu()
{
}




template<typename T, typename I>
void TotalFETIExplicitGeneralScCpu<T,I>::info()
{
    eslog::info(" = EXPLICIT TOTAL FETI OPERATOR USING TRIANGULAR SC                                   = \n");
    eslog::info(" =   EXTERNAL SPARSE SOLVER               %50s = \n", DirectSparseSolver<T>::name());

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
void TotalFETIExplicitGeneralScCpu<T,I>::set(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralScCpu::set");

    n_domains = feti.K.size();

    domain_data.resize(n_domains);

    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        data.n_dofs_domain = feti.B1[di].ncols;
        data.n_dofs_interface = feti.B1[di].nrows;
    }

    replace_unset_configs<T,I>(cfg);

    Fs_allocated.resize((n_domains - 1) / 2 + 1);
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
            MatrixDenseData_new<T> & F_allocd = Fs_allocated[allocated_F_index];
            per_domain_stuff & data_bigger = domain_data[di_bigger];
            per_domain_stuff & data_di = domain_data[di];
            if(di == di_bigger) {
                F_allocd.set(data_bigger.n_dofs_interface + 1, data_bigger.n_dofs_interface, cfg.order_F, AllocatorCPU_new::get_singleton());
                F_allocd.alloc();
                data_di.F = F_allocd.get_submatrix_view(1, data_di.n_dofs_interface + 1, 0, data_di.n_dofs_interface);
                data_di.F.prop.uplo = 'L';
            }
            else {
                data_di.F = F_allocd.get_submatrix_view(0, data_di.n_dofs_interface , 0, data_di.n_dofs_interface);
                data_di.F.prop.uplo = 'U';
            }
            data_di.F.prop.symm = MatrixSymmetry_new::hermitian;
            if(data_di.F.order == 'R') {
                data_di.F_old = MatrixDenseView_new<T>::template to_old<I,cpu_allocator>(data_di.F);
                data_di.F_uplo_in_rowmajor = data_di.F.prop.uplo;
            }
            if(data_di.F.order == 'C') {
                MatrixDenseView_new<T> F_reordered = data_di.F.get_transposed_reordered_view();
                data_di.F_old = MatrixDenseView_new<T>::template to_old<I,cpu_allocator>(F_reordered);
                data_di.F_uplo_in_rowmajor = change_uplo(data_di.F.prop.uplo);
            }
        }
    }

    stacktimer::push("TotalFETIExplicitGeneralScCpu::set preprocess");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];

        stacktimer::info("TotalFETIExplicitGeneralScCpu::set preprocess subdomain %zu", di);

        math::combine(data.Kreg_old, feti.K[di], feti.RegMat[di]);
        if constexpr(utils::is_real<T>())    data.Kreg_old.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
        if constexpr(utils::is_complex<T>()) data.Kreg_old.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE;
        data.Kreg_old.shape = feti.K[di].shape;

        data.Bt = MatrixCsxView_new<T,I>::from_old(feti.B1[di]).get_transposed_reordered_view();
        data.Kreg = MatrixCsxView_new<T,I>::from_old(data.Kreg_old);
        if constexpr(utils::is_real<T>()) if(is_symmetric<T>(data.Kreg.prop.symm)) data.Kreg.prop.symm = MatrixSymmetry_new::hermitian;

        data.op_sc = math::operations::sc_csx_dny<T,I>::make();
        data.op_sc->set_coefficients(-1);
        data.op_sc->set_matrix(&data.Kreg, &data.Bt, nullptr, nullptr);
        data.op_sc->set_sc(&data.F);
        data.op_sc->set_need_solve_A11(true);
        data.op_sc->preprocess();

        data.x.resize(data.n_dofs_interface);
        data.y.resize(data.n_dofs_interface);
    }
    if(!cfg.inner_timers) stacktimer::enable();
    stacktimer::pop();

    // clean up the mess from buggy openmp in clang
    utils::run_dummy_parallel_region();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIExplicitGeneralScCpu<T,I>::update(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralScCpu::update");

    stacktimer::push("update_mainloop");
    if(cfg.mainloop_update_split == 'C') {
        stacktimer::push("update_mainloop_combined");
        if(!cfg.inner_timers) stacktimer::disable();
        #pragma omp parallel for schedule(static,1) if(cfg.parallel_update)
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_stuff & data = domain_data[di];
            math::sumCombined(data.Kreg_old, T{1.0}, feti.K[di], feti.RegMat[di]);
            data.op_sc->perform();
        };
        if(!cfg.inner_timers) stacktimer::enable();
        stacktimer::pop();

        // clean up the mess from buggy openmp in clang
        utils::run_dummy_parallel_region();
    }
    if(cfg.mainloop_update_split == 'S') {
        stacktimer::push("update_mainloop_separate_1");
        if(!cfg.inner_timers) stacktimer::disable();
        #pragma omp parallel for schedule(static,1) if(cfg.parallel_update)
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_stuff & data = domain_data[di];
            math::sumCombined(data.Kreg_old, T{1.0}, feti.K[di], feti.RegMat[di]);
            data.op_sc->perform_1();
        }
        if(!cfg.inner_timers) stacktimer::enable();
        stacktimer::pop();

        // clean up the mess from buggy openmp in clang
        utils::run_dummy_parallel_region();

        stacktimer::push("update_mainloop_separate_2");
        if(!cfg.inner_timers) stacktimer::disable();
        #pragma omp parallel for schedule(static,1) if(cfg.parallel_update)
        for(size_t di = 0; di < n_domains; di++) {
            domain_data[di].op_sc->perform_2();
        }
        if(!cfg.inner_timers) stacktimer::enable();
        stacktimer::pop();
    }
    stacktimer::pop();

    stacktimer::push("update_compute_vector_d");
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
    stacktimer::pop();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIExplicitGeneralScCpu<T,I>::_apply(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralScCpu::apply");

    std::fill_n(y_cluster.vals, y_cluster.size, T{0});

    #pragma omp parallel for schedule(static,1) if(cfg.parallel_apply)
    for(size_t di = 0; di < n_domains; di++) {
        auto & data = domain_data[di];

        std::vector<I> & D2C = feti.D2C[di];

        for(size_t i = 0; i < data.n_dofs_interface; i++) {
            data.x.vals[i] = x_cluster.vals[D2C[i]];
        }

        math::blas::apply_hermitian(data.y, T{1}, data.F_old, data.F_uplo_in_rowmajor, T{0}, data.x);

        for(size_t i = 0; i < data.n_dofs_interface; i++) {
            T val = data.y.vals[i];
            T & dst = y_cluster.vals[D2C[i]];
            if constexpr(utils::is_complex<T>()) {
                #pragma omp atomic
                utils::real_ref(dst) += utils::real_ref(val);
                #pragma omp atomic
                utils::imag_ref(dst) += utils::imag_ref(val);
            }
            else {
                #pragma omp atomic
                dst += val;
            }
        }
    }

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIExplicitGeneralScCpu<T,I>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    _apply(x, y);
    y.synchronize();
}



template<typename T, typename I>
void TotalFETIExplicitGeneralScCpu<T,I>::apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y)
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
void TotalFETIExplicitGeneralScCpu<T,I>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
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
}



template<typename T, typename I>
void TotalFETIExplicitGeneralScCpu<T,I>::print(const step::Step &step)
{
    eslog::error("TotalFETIExplicitGeneralScCpu::print not implemented");
}



#define INSTANTIATE_T_I(T,I) \
template class TotalFETIExplicitGeneralScCpu<T,I>;

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

