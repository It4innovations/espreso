
#include "feti/dualoperator/totalfeti.implicit.generalsparsesolver.cpu.h"

#include "math/primitives_new/allocator_new.h"
#include "feti/common/applyB.h"
#include "math/math.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/minmaxavg.h"
#include "basis/utilities/stacktimer.h"
#include "math/operations/gemv_csx.h"
#include "math/operations/gemm_csx_dny_dnz.h"
#include "math/operations/fill_dnx.h"
#include "math/operations/submatrix_dnx_dnx_noncontig.h"
#include "math/operations/supermatrix_dnx_dnx_noncontig.h"



namespace espreso {



template<typename T, typename I>
TotalFETIImplicitGeneralSparseSolverCpu<T,I>::TotalFETIImplicitGeneralSparseSolverCpu(FETI<T> &feti) : DualOperator<T>(feti)
{
    setup_config(cfg, feti.configuration);
}



template<typename T, typename I>
TotalFETIImplicitGeneralSparseSolverCpu<T,I>::~TotalFETIImplicitGeneralSparseSolverCpu()
{
}



template<typename T, typename I>
void TotalFETIImplicitGeneralSparseSolverCpu<T,I>::info()
{
    eslog::info(" = IMPLICIT TOTAL FETI OPERATOR USING GENERAL SPARSE SOLVER ON CPU                           = \n");

    if(cfg.print_config) {
        auto bool_to_string = [](bool val){ return val ? "TRUE" : "FALSE";};
        auto solver_impl_to_string = [](solver_impl_t solver_impl){ switch(solver_impl) { case solver_impl_t::autoselect: return "autoselect"; case solver_impl_t::mklpardiso: return "mklpardiso"; case solver_impl_t::suitesparse: return "suitesparse"; case solver_impl_t::mumps: return "mumps"; case solver_impl_t::strumpack: return "strumpack"; case solver_impl_t::pastix: return "pastix"; case solver_impl_t::superlu_dist: return "superlu_dist"; default: return "UNDEFINED"; } };
        auto solver_impl_to_string_actual = [&](solver_impl_t solver_impl){
            if(feti.K.size() == 0) return "UNDEFINED";
            MatrixBase_new::matrix_properties prop;
            prop.symm = get_new_matrix_symmetry(feti.K[0].type);
            prop.dfnt = get_new_matrix_definitness(feti.K[0].type);
            return math::operations::solver_csx<T,I>::make(solver_impl, &prop, false, true)->get_name();
        };

        eslog::info(" =   %-50s       %+30s = \n", "parallel_set", bool_to_string(cfg.parallel_set));
        eslog::info(" =   %-50s       %+30s = \n", "parallel_update", bool_to_string(cfg.parallel_update));
        eslog::info(" =   %-50s       %+30s = \n", "parallel_apply", bool_to_string(cfg.parallel_apply));
        eslog::info(" =   %-50s       %+30s = \n", "inner_timers", bool_to_string(cfg.inner_timers));
        eslog::info(" =   %-50s       %+30s = \n", "outer_timers", bool_to_string(cfg.outer_timers));
        eslog::info(" =   %-50s       %+30s = \n", "solver_impl", solver_impl_to_string(cfg.solver_impl));
        eslog::info(" =   %-50s       %+30s = \n", "solver_impl_actual", solver_impl_to_string_actual(cfg.solver_impl));
    }

    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_domain; }).to_string("  Domain volume [dofs]").c_str());
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_interface; }).to_string("  Domain interface [dofs]").c_str());
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}



template<typename T, typename I>
void TotalFETIImplicitGeneralSparseSolverCpu<T,I>::set(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIImplicitGeneralSparseSolverCpu::set");

    n_domains = feti.K.size();

    domain_data.resize(n_domains);

    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & my = domain_data[di];
        my.n_dofs_domain = feti.B1[di].ncols;
        my.n_dofs_interface = feti.B1[di].nrows;
    }

    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & my = domain_data[di];

        stacktimer::info("TotalFETIExplicitGeneralScCpu::set preprocess subdomain %zu", di);

        math::combine(my.Kreg_old, feti.K[di], feti.RegMat[di]);
        if constexpr(utils::is_real<T>())    my.Kreg_old.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
        if constexpr(utils::is_complex<T>()) my.Kreg_old.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE;
        my.Kreg_old.shape = feti.K[di].shape;

        my.Kreg = MatrixCsxView_new<T,I>::from_old(my.Kreg_old);

        my.op_solver = math::operations::solver_csx<T,I>::make(cfg.solver_impl, &my.Kreg.prop, false, true);
        my.op_solver->set_matrix_A(&my.Kreg);
        my.op_solver->set_needs(false, true);
        my.op_solver->factorize_symbolic();
    }
    if(!cfg.inner_timers) stacktimer::enable();

    // clean up the mess from buggy openmp in clang
    utils::run_dummy_parallel_region();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIImplicitGeneralSparseSolverCpu<T,I>::update(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIImplicitGeneralSparseSolverCpu::update");

    stacktimer::push("update_factorize_numeric");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_update)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & my = domain_data[di];

        math::sumCombined(my.Kreg_old, T{1.0}, feti.K[di], feti.RegMat[di]);
        my.op_solver->factorize_numeric();
    };
    if(!cfg.inner_timers) stacktimer::enable();
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
            domain_data[di].op_solver->solve(f_new, Kplus_f_new);
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
void TotalFETIImplicitGeneralSparseSolverCpu<T,I>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIImplicitGeneralSparseSolverCpu::apply (vector)");

    stacktimer::push("cpu_implicit_apply_vector_compute");

    std::fill_n(y.vals, y.size, T{0});

    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_apply)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & my = domain_data[di];
        int * D2C = feti.D2C[di].data();

        VectorDenseData_new<T> z;
        z.set(my.n_dofs_interface, AllocatorCPU_new::get_singleton());
        z.alloc();
        
        VectorDenseData_new<T> w;
        w.set(my.n_dofs_domain, AllocatorCPU_new::get_singleton());
        w.alloc();

        for(size_t i = 0; i < my.n_dofs_interface; i++) {
            z.vals[i] = x.vals[D2C[i]];
        }

        MatrixCsxView_new<T,I> B = MatrixCsxView_new<T,I>::from_old(feti.B1[di]);
        MatrixCsxView_new<T,I> Bt = B.get_transposed_reordered_view();

        math::operations::gemv_csx<T,I>::do_all(&Bt, &z, &w, T{1}, T{0});

        my.op_solver->solve(w, w);

        math::operations::gemv_csx<T,I>::do_all(&B, &w, &z, T{1}, T{0});

        for(size_t i = 0; i < my.n_dofs_interface; i++) {
            utils::atomic_add(y.vals[D2C[i]], z.vals[i]);
        }
    }
    if(!cfg.inner_timers) stacktimer::enable();

    stacktimer::pop();

    y.synchronize();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIImplicitGeneralSparseSolverCpu<T,I>::apply(const Matrix_Dual<T> &X_cluster_old, Matrix_Dual<T> &Y_cluster_old)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIImplicitGeneralSparseSolverCpu::apply (matrix)");

    MatrixDenseView_new<T> X_cluster = MatrixDenseView_new<T>::from_old(X_cluster_old).get_transposed_reordered_view();
    MatrixDenseView_new<T> Y_cluster = MatrixDenseView_new<T>::from_old(Y_cluster_old).get_transposed_reordered_view();

    size_t n_dofs_cluster_interface = feti.lambdas.size;
    if(X_cluster.nrows != Y_cluster.nrows || X_cluster.ncols != Y_cluster.ncols) eslog::error("size of X_cluster and Y_cluster does not match\n");
    if(X_cluster.order != Y_cluster.order) eslog::error("orders do not match\n");
    if(X_cluster.nrows != n_dofs_cluster_interface) eslog::error("incompatible cluster matrix size\n");
    char order = X_cluster.order;
    size_t width = X_cluster.ncols;

    stacktimer::push("cpu_implicit_apply_matrix_compute");

    math::operations::fill_dnx<T>::do_all(&Y_cluster, T{0});

    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_apply)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & my = domain_data[di];

        VectorDenseView_new<I> my_D2C;
        my_D2C.set_view(feti.D2C[di].size(), feti.D2C[di].data(), AllocatorDummy_new::get_singleton(true, false));

        MatrixDenseData_new<T> Z;
        Z.set(my.n_dofs_interface, width, order, AllocatorCPU_new::get_singleton());
        Z.alloc();

        MatrixDenseData_new<T> W;
        W.set(my.n_dofs_domain, width, order, AllocatorCPU_new::get_singleton());
        W.alloc();

        math::operations::submatrix_dnx_dnx_noncontig<T,I>::do_all(&X_cluster, &Z, &my_D2C, nullptr);

        MatrixCsxView_new<T,I> B = MatrixCsxView_new<T,I>::from_old(feti.B1[di]);
        MatrixCsxView_new<T,I> Bt = B.get_transposed_reordered_view();

        math::operations::gemm_csx_dny_dnz<T,I>::do_all(&Bt, &Z, &W, T{1}, T{0});

        my.op_solver->solve(W, W);

        math::operations::gemm_csx_dny_dnz<T,I>::do_all(&B, &W, &Z, T{1}, T{0});

        math::operations::supermatrix_dnx_dnx_noncontig<T,I>::do_all(&Z, &Y_cluster, &my_D2C, nullptr, math::operations::supermatrix_dnx_dnx_noncontig<T,I>::mode::accumulate_atomic);
    }
    if(!cfg.inner_timers) stacktimer::enable();

    stacktimer::pop();
    
    Y_cluster_old.synchronize();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIImplicitGeneralSparseSolverCpu<T,I>::apply(const Matrix_Dual<T> &X_cluster_old, Matrix_Dual<T> &Y_cluster_old, const std::vector<int> &filter)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIImplicitGeneralSparseSolverCpu::apply (matrix,filter)");

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

    size_t n_dofs_cluster_interface = feti.lambdas.size;
    if(X_cluster.nrows != Y_cluster.nrows || X_cluster.ncols != Y_cluster.ncols) eslog::error("size of X_cluster and Y_cluster does not match\n");
    if(X_cluster.order != Y_cluster.order) eslog::error("orders do not match\n");
    if(X_cluster.nrows != n_dofs_cluster_interface) eslog::error("incompatible cluster matrix size\n");
    char order = X_cluster.order;
    size_t width = filter_map.size;

    stacktimer::push("cpu_implicit_apply_matrix_compute");

    math::operations::fill_dnx<T>::do_all(&Y_cluster, T{0});

    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_apply)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & my = domain_data[di];

        VectorDenseView_new<I> my_D2C;
        my_D2C.set_view(feti.D2C[di].size(), feti.D2C[di].data(), AllocatorDummy_new::get_singleton(true, false));

        MatrixDenseData_new<T> Z;
        Z.set(my.n_dofs_interface, width, order, AllocatorCPU_new::get_singleton());
        Z.alloc();

        MatrixDenseData_new<T> W;
        W.set(my.n_dofs_domain, width, order, AllocatorCPU_new::get_singleton());
        W.alloc();

        math::operations::submatrix_dnx_dnx_noncontig<T,I>::do_all(&X_cluster, &Z, &my_D2C, &filter_map);

        MatrixCsxView_new<T,I> B = MatrixCsxView_new<T,I>::from_old(feti.B1[di]);
        MatrixCsxView_new<T,I> Bt = B.get_transposed_reordered_view();

        math::operations::gemm_csx_dny_dnz<T,I>::do_all(&Bt, &Z, &W, T{1}, T{0});

        my.op_solver->solve(W, W);

        math::operations::gemm_csx_dny_dnz<T,I>::do_all(&B, &W, &Z, T{1}, T{0});

        math::operations::supermatrix_dnx_dnx_noncontig<T,I>::do_all(&Z, &Y_cluster, &my_D2C, &filter_map, math::operations::supermatrix_dnx_dnx_noncontig<T,I>::mode::accumulate_atomic);
    }
    if(!cfg.inner_timers) stacktimer::enable();

    stacktimer::pop();
    
    Y_cluster_old.synchronize();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIImplicitGeneralSparseSolverCpu<T,I>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIImplicitGeneralSparseSolverCpu::toPrimal");
    if(!cfg.inner_timers) stacktimer::disable();

    #pragma omp parallel for schedule(static,1)
    for (size_t di = 0; di < n_domains; ++di) {
        Vector_Dense<T,I> z;
        z.resize(y[di]);
        applyBt(feti, di, x, z, T{-1});
        math::add(z, T{1}, feti.f[di]);
        VectorDenseView_new<T> z_new = VectorDenseView_new<T>::from_old(z);
        VectorDenseView_new<T> y_new = VectorDenseView_new<T>::from_old(y[di]);
        domain_data[di].op_solver->solve(z_new, y_new);
    }

    if(!cfg.inner_timers) stacktimer::enable();
    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIImplicitGeneralSparseSolverCpu<T,I>::setup_config(config & cfg, const FETIConfiguration & feti_ecf_config)
{
    // defaults are set in config definition
    // if ecf value is auto, cfg value is not changed

    using ecf_config = DualopTotalfetiImplicitGeneralSparseSolverCpuConfig;
    const ecf_config & ecf = feti_ecf_config.dualop_totalfeti_implicit_generalsparsesolver_cpu_config;

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

    switch(ecf.sparse_solver_impl) {
        case ecf_config::SPARSE_SOLVER_IMPL::AUTO: break;
        case ecf_config::SPARSE_SOLVER_IMPL::MKLPARDISO:   cfg.solver_impl = solver_impl_t::mklpardiso;   break;
        case ecf_config::SPARSE_SOLVER_IMPL::SUITESPARSE:  cfg.solver_impl = solver_impl_t::suitesparse;  break;
        case ecf_config::SPARSE_SOLVER_IMPL::MUMPS:        cfg.solver_impl = solver_impl_t::mumps;        break;
        case ecf_config::SPARSE_SOLVER_IMPL::STRUMPACK:    cfg.solver_impl = solver_impl_t::strumpack;    break;
        case ecf_config::SPARSE_SOLVER_IMPL::PASTIX:       cfg.solver_impl = solver_impl_t::pastix;       break;
        case ecf_config::SPARSE_SOLVER_IMPL::SUPERLU_DIST: cfg.solver_impl = solver_impl_t::superlu_dist; break;
    }
}



#define INSTANTIATE_T_I(T,I) \
template class TotalFETIImplicitGeneralSparseSolverCpu<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

        #define INSTANTIATE \
        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        INSTANTIATE_T(std::complex<double>)

            INSTANTIATE

    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
