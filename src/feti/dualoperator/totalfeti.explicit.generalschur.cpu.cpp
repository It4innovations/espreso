
#include "feti/dualoperator/totalfeti.explicit.generalschur.cpu.h"

#include "math/primitives_new/allocator_new.h"
#include "feti/common/applyB.h"
#include "basis/utilities/minmaxavg.h"
#include "basis/utilities/stacktimer.h"
#include "esinfo/meshinfo.h"

#include <algorithm>



namespace espreso {



template<typename T, typename I>
TotalFETIExplicitGeneralSchurCpu<T,I>::TotalFETIExplicitGeneralSchurCpu(FETI<T> &feti) : DualOperator<T>(feti)
{
    setup_config(cfg, feti.configuration);
}



template<typename T, typename I>
TotalFETIExplicitGeneralSchurCpu<T,I>::~TotalFETIExplicitGeneralSchurCpu()
{
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurCpu<T,I>::info()
{
    eslog::info(" = EXPLICIT TOTAL FETI OPERATOR USING GENERAL SCHUR ON CPU                                   = \n");

    if(cfg.print_config) {
        auto order_to_string = [](char order){ switch(order){ case 'R': return "ROW_MAJOR"; case 'C': return "COL_MAJOR"; default: return "UNDEFINED"; }};
        auto bool_to_string = [](bool val){ return val ? "TRUE" : "FALSE";};
        auto loop_split_to_string = [](char val){ switch(val){ case 'C': return "COMBINED"; case 'S': return "SEPARATE"; default: return "UNDEFINED"; }};
        auto schur_impl_to_string = [](schur_impl_t schur_impl){ switch(schur_impl) { case schur_impl_t::autoselect: return "autoselect"; case schur_impl_t::triangular: return "triangular"; case schur_impl_t::mklpardiso: return "mklpardiso"; case schur_impl_t::sparse_solver: return "sparse_solver"; case schur_impl_t::mumps: return "mumps"; case schur_impl_t::pastix: return "pastix"; default: return "UNDEFINED"; }};
        auto schur_impl_to_string_actual = [](schur_impl_t schur_impl){ return math::operations::schur_csx_dny<T,I>::make(schur_impl)->get_name(); };

        eslog::info(" =   %-50s       %+30s = \n", "parallel_set", bool_to_string(cfg.parallel_set));
        eslog::info(" =   %-50s       %+30s = \n", "parallel_update", bool_to_string(cfg.parallel_update));
        eslog::info(" =   %-50s       %+30s = \n", "parallel_apply", bool_to_string(cfg.parallel_apply));
        eslog::info(" =   %-50s       %+30s = \n", "mainloop_update_split", loop_split_to_string(cfg.mainloop_update_split));
        eslog::info(" =   %-50s       %+30s = \n", "inner_timers", bool_to_string(cfg.inner_timers));
        eslog::info(" =   %-50s       %+30s = \n", "outer_timers", bool_to_string(cfg.outer_timers));
        eslog::info(" =   %-50s       %+30s = \n", "order_F", order_to_string(cfg.order_F));
        eslog::info(" =   %-50s       %+30s = \n", "schur implementation", schur_impl_to_string(cfg.schur_impl));
        eslog::info(" =   %-50s       %+30s = \n", "schur implementation actual", schur_impl_to_string_actual(cfg.schur_impl));
    }

    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_domain; }).to_string("  Domain volume [dofs]").c_str());
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_interface; }).to_string("  Domain surface [dofs]").c_str());
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurCpu<T,I>::set(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralSchurCpu::set");

    n_domains = feti.K.size();

    domain_data.resize(n_domains);

    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        data.n_dofs_domain = feti.B1[di].ncols;
        data.n_dofs_interface = feti.B1[di].nrows;
    }

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
                F_allocd.set(data_bigger.n_dofs_interface + (cfg.order_F == 'R'), data_bigger.n_dofs_interface + (cfg.order_F == 'C'), cfg.order_F, AllocatorCPU_new::get_singleton());
                F_allocd.alloc();
                data_di.F.set_view(data_di.n_dofs_interface, data_di.n_dofs_interface, F_allocd.ld, F_allocd.order, F_allocd.vals + F_allocd.ld, F_allocd.ator);
                data_di.F.prop.uplo = (F_allocd.order == 'R') ? 'L' : 'U';
            }
            else {
                data_di.F = F_allocd.get_submatrix_view(0, data_di.n_dofs_interface, 0, data_di.n_dofs_interface);
                data_di.F.prop.uplo = (F_allocd.order == 'R') ? 'U' : 'L';
            }
            data_di.F.prop.symm = MatrixSymmetry_new::hermitian;
            if(data_di.F.order == 'R') {
                data_di.F_old = MatrixDenseView_new<T>::template to_old<I,cpu_allocator>(data_di.F);
            }
            if(data_di.F.order == 'C') {
                MatrixDenseView_new<T> F_reordered = data_di.F.get_transposed_reordered_view();
                data_di.F_old = MatrixDenseView_new<T>::template to_old<I,cpu_allocator>(F_reordered);
            }
        }
    }

    {
        std::vector<MatrixDenseView_new<T>*> Fs_vector(n_domains);
        for(size_t di = 0; di < n_domains; di++) {
            Fs_vector[di] = &domain_data[di].F;
        }

        applicator.set_handles(&feti.main_q, &feti.queues, &feti.handles_dense);
        applicator.set_dimensions(feti);
        applicator.set_vector_memory('C');
        applicator.set_D2C_map(&feti.D2C);
        applicator.set_Fs(Fs_vector);
        applicator.set_apply_target('C');
        applicator.setup();
        applicator.preprocess();
    }

    stacktimer::push("TotalFETIExplicitGeneralSchurCpu::set preprocess");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];

        stacktimer::info("TotalFETIExplicitGeneralSchurCpu::set preprocess subdomain %zu", di);

        math::combine(data.Kreg_old, feti.K[di], feti.RegMat[di]);
        if constexpr(utils::is_real<T>())    data.Kreg_old.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
        if constexpr(utils::is_complex<T>()) data.Kreg_old.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE;
        data.Kreg_old.shape = feti.K[di].shape;

        data.Bt = MatrixCsxView_new<T,I>::from_old(feti.B1[di]).get_transposed_reordered_view();
        data.Kreg = MatrixCsxView_new<T,I>::from_old(data.Kreg_old);
        if constexpr(utils::is_real<T>()) if(is_symmetric<T>(data.Kreg.prop.symm)) data.Kreg.prop.symm = MatrixSymmetry_new::hermitian;

        data.op_sc = math::operations::schur_csx_dny<T,I>::make(cfg.schur_impl);
        data.op_sc->set_coefficients(-1);
        data.op_sc->set_matrix(&data.Kreg, &data.Bt, nullptr, nullptr);
        data.op_sc->set_sc(&data.F);
        data.op_sc->set_need_solve_A11(true);
        data.op_sc->preprocess();
    }
    if(!cfg.inner_timers) stacktimer::enable();
    stacktimer::pop();

    // clean up the mess from buggy openmp in clang
    utils::run_dummy_parallel_region();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurCpu<T,I>::update(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralSchurCpu::update");

    stacktimer::push("update_mainloop");
    if(cfg.mainloop_update_split == 'C') {
        stacktimer::push("update_mainloop_combined");
        if(!cfg.inner_timers) stacktimer::disable();
        #pragma omp parallel for schedule(static,1) if(cfg.parallel_update)
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_stuff & data = domain_data[di];
            math::sumCombined(data.Kreg_old, T{1.0}, feti.K[di], feti.RegMat[di]);
            data.op_sc->perform();
            applicator.update_F(di);
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
            applicator.update_F(di);
        }
        if(!cfg.inner_timers) stacktimer::enable();
        stacktimer::pop();
    }
    stacktimer::pop();

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

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurCpu<T,I>::_apply(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    VectorDenseView_new<T> x_cluster_new = VectorDenseView_new<T>::from_old(x_cluster);
    VectorDenseView_new<T> y_cluster_new = VectorDenseView_new<T>::from_old(y_cluster);

    applicator.apply(x_cluster_new, y_cluster_new);
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurCpu<T,I>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralSchurCpu::apply");

    _apply(x, y);
    y.synchronize();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurCpu<T,I>::apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralSchurCpu::apply");

    Vector_Dual<T> _x, _y;
    _x.size = _y.size = x.ncols;
    for (int r = 0; r < x.nrows; ++r) {
        _x.vals = x.vals + x.ncols * r;
        _y.vals = y.vals + y.ncols * r;
        _apply(_x, _y);
    }
    y.synchronize();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurCpu<T,I>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIExplicitGeneralSchurCpu::toPrimal");
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
void TotalFETIExplicitGeneralSchurCpu<T,I>::print(const step::Step &step)
{
    eslog::error("TotalFETIExplicitGeneralSchurCpu::print not implemented");
}



template<typename T, typename I>
void TotalFETIExplicitGeneralSchurCpu<T,I>::setup_config(config & cfg, const FETIConfiguration & feti_ecf_config)
{
    // defaults are set in config definition
    // if ecf value is auto, cfg value is not changed

    using ecf_config = DualopTotalfetiExplicitGeneralSchurCpuConfig;
    const ecf_config & ecf = feti_ecf_config.dualop_totalfeti_explicit_generalschur_cpu_config;

    if(ecf.parallel_set == ecf_config::AUTOBOOL::TRUE)  cfg.parallel_set = true;
    if(ecf.parallel_set == ecf_config::AUTOBOOL::FALSE) cfg.parallel_set = false;

    if(ecf.parallel_update == ecf_config::AUTOBOOL::TRUE)  cfg.parallel_update = true;
    if(ecf.parallel_update == ecf_config::AUTOBOOL::FALSE) cfg.parallel_update = false;

    if(ecf.parallel_apply == ecf_config::AUTOBOOL::TRUE)  cfg.parallel_apply = true;
    if(ecf.parallel_apply == ecf_config::AUTOBOOL::FALSE) cfg.parallel_apply = false;

    if(ecf.mainloop_update_split == ecf_config::MAINLOOP_UPDATE_SPLIT::COMBINED) cfg.mainloop_update_split = 'C';
    if(ecf.mainloop_update_split == ecf_config::MAINLOOP_UPDATE_SPLIT::SEPARATE) cfg.mainloop_update_split = 'S';

    if(ecf.timers_outer == ecf_config::AUTOBOOL::TRUE)  cfg.outer_timers = true;
    if(ecf.timers_outer == ecf_config::AUTOBOOL::FALSE) cfg.outer_timers = false;

    if(ecf.timers_inner == ecf_config::AUTOBOOL::TRUE)  cfg.inner_timers = true;
    if(ecf.timers_inner == ecf_config::AUTOBOOL::FALSE) cfg.inner_timers = false;

    if(ecf.print_config == ecf_config::AUTOBOOL::TRUE)  cfg.print_config = true;
    if(ecf.print_config == ecf_config::AUTOBOOL::FALSE) cfg.print_config = false;

    if(ecf.order_F == ecf_config::MATRIX_ORDER::ROW_MAJOR) cfg.order_F = 'R';
    if(ecf.order_F == ecf_config::MATRIX_ORDER::COL_MAJOR) cfg.order_F = 'C';

    if(ecf.schur_impl == ecf_config::SCHUR_IMPL::TRIANGULAR)    cfg.schur_impl = schur_impl_t::triangular;
    if(ecf.schur_impl == ecf_config::SCHUR_IMPL::MKLPARDISO)    cfg.schur_impl = schur_impl_t::mklpardiso;
    if(ecf.schur_impl == ecf_config::SCHUR_IMPL::SPARSE_SOLVER) cfg.schur_impl = schur_impl_t::sparse_solver;
    if(ecf.schur_impl == ecf_config::SCHUR_IMPL::MUMPS)         cfg.schur_impl = schur_impl_t::mumps;
    if(ecf.schur_impl == ecf_config::SCHUR_IMPL::PASTIX)        cfg.schur_impl = schur_impl_t::pastix;
}



#define INSTANTIATE_T_I(T,I) \
template class TotalFETIExplicitGeneralSchurCpu<T,I>;

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

