
#include "feti/dualoperator/hybridfeti.implicit.generalsparsesolver.cpu.h"

#include "math/primitives_new/allocator_new.h"
#include "feti/common/applyB.h"
#include "math/math.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/clusterstore.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/sysutils.h"
#include "basis/utilities/minmaxavg.h"
#include "basis/utilities/stacktimer.h"
#include "math/operations/gemv_csx.h"
#include "math/operations/gemm_csx_dny_dnz.h"
#include "math/operations/fill_dnx.h"
#include "math/operations/submatrix_dnx_dnx_noncontig.h"
#include "math/operations/supermatrix_dnx_dnx_noncontig.h"



namespace espreso {



template<typename T, typename I>
HybridFETIImplicitGeneralSparseSolverCpu<T,I>::HybridFETIImplicitGeneralSparseSolverCpu(FETI<T> &feti) : DualOperator<T>(feti)
{
    setup_config(cfg, feti.configuration);
}



template<typename T, typename I>
HybridFETIImplicitGeneralSparseSolverCpu<T,I>::~HybridFETIImplicitGeneralSparseSolverCpu()
{
}



template<typename T, typename I>
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::info()
{
    eslog::info(" = IMPLICIT HYBRID FETI OPERATOR USING GENERAL SPARSE SOLVER ON CPU                          = \n");

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
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::set(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIImplicitGeneralSparseSolverCpu::set");

    n_domains = feti.K.size();

    Btx.resize(feti.K.size());
    KplusBtx.resize(feti.K.size());
    dKB0.resize(feti.K.size());

    domain_data.resize(n_domains);

    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & my = domain_data[di];
        my.n_dofs_domain = feti.B1[di].ncols;
        my.n_dofs_interface = feti.B1[di].nrows;
    }

    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & my = domain_data[di];

        stacktimer::info("preprocess subdomain %zu", di);

        math::combine(my.Kreg_old, feti.K[di], feti.RegMat[di]);
        if constexpr(utils::is_real<T>())    my.Kreg_old.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
        if constexpr(utils::is_complex<T>()) my.Kreg_old.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE;
        my.Kreg_old.shape = feti.K[di].shape;

        my.Kreg = MatrixCsxView_new<T,I>::from_old(my.Kreg_old);

        my.op_solver = math::operations::solver_csx<T,I>::make(cfg.solver_impl, &my.Kreg.prop, false, true);
        my.op_solver->set_matrix_A(&my.Kreg);
        my.op_solver->set_needs(false, true);
        my.op_solver->factorize_symbolic();

        Btx[di].resize(feti.K[di].nrows);
        KplusBtx[di].resize(feti.K[di].nrows);
        math::set(Btx[di], T{0});
    }
    if(!cfg.inner_timers) stacktimer::enable();

    // clean up the mess from buggy openmp in clang
    utils::run_dummy_parallel_region();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::update(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIImplicitGeneralSparseSolverCpu::update");

    stacktimer::push("update_factorize_numeric");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & my = domain_data[di];

        math::sumCombined(my.Kreg_old, T{1.0}, feti.K[di], feti.RegMat[di]);
        my.op_solver->factorize_numeric();
    };
    if(!cfg.inner_timers) stacktimer::enable();
    stacktimer::pop();

    stacktimer::push("update_second_part");
    {
        stacktimer::push("hybrid_feti_things");
        if(!cfg.inner_timers) stacktimer::disable();
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
        if(!cfg.inner_timers) stacktimer::enable();
        stacktimer::pop();

        stacktimer::push("update_compute_vector_d");
        if(!cfg.inner_timers) stacktimer::disable();
        {
            if (feti.updated.B) {
                d.resize();
            }
            _applyK(feti.f, KplusBtx, true);
            applyB(feti, KplusBtx, d);
            d.synchronize();
            math::add(d, T{-1}, feti.c);
            eslog::checkpointln("FETI: COMPUTE DUAL RHS [d]");
            if (info::ecf->output.print_matrices) {
                eslog::storedata(" STORE: feti/dualop/{d}\n");
                math::store(d, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "d").c_str());
            }
        }
        if(!cfg.inner_timers) stacktimer::enable();
        stacktimer::pop();
    }
    stacktimer::pop();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIImplicitGeneralSparseSolverCpu::apply (vector)");

    stacktimer::push("cpu_implicit_apply_vector_compute");

    std::fill_n(y.vals, y.size, T{0});

    _apply_hfeti_stuff(x, y);

    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1)
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

    stacktimer::push("dual_synchronize");
    y.synchronize();
    stacktimer::pop();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::apply(const Matrix_Dual<T> &X_cluster_old, Matrix_Dual<T> &Y_cluster_old)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIImplicitGeneralSparseSolverCpu::apply (matrix)");

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

    _apply_hfeti_stuff(X_cluster_old, Y_cluster_old);

    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1)
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
    
    stacktimer::push("dual_synchronize");
    Y_cluster_old.synchronize();
    stacktimer::pop();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::apply(const Matrix_Dual<T> &X_cluster_old, Matrix_Dual<T> &Y_cluster_old, const std::vector<int> &filter)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIImplicitGeneralSparseSolverCpu::apply (matrix,filter)");

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

    _apply_hfeti_stuff(X_cluster_old, Y_cluster_old, filter);

    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1)
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
    
    stacktimer::push("dual_synchronize");
    Y_cluster_old.synchronize();
    stacktimer::pop();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::apply(const Matrix_Dual<T> &X_cluster_old, Matrix_Dual<T> &Y_cluster_old, const std::vector<std::vector<int>> &filter)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIImplicitGeneralSparseSolverCpu::apply (matrix,filter2)");

    MatrixDenseView_new<T> X_cluster = MatrixDenseView_new<T>::from_old(X_cluster_old).get_transposed_reordered_view();
    MatrixDenseView_new<T> Y_cluster = MatrixDenseView_new<T>::from_old(Y_cluster_old).get_transposed_reordered_view();

    // orig filter: i-th rhs is handled by onyl a subset of domains
    // transposed: i-th domain handles only a subset of rhs
    std::vector<std::vector<I>> filter_transposed(n_domains);
    for(size_t idx_rhs = 0; idx_rhs < filter.size(); idx_rhs++) {
        for(int di : filter[idx_rhs]) {
            filter_transposed[di].push_back(idx_rhs);
        }
    }

    size_t n_dofs_cluster_interface = feti.lambdas.size;
    if(X_cluster.nrows != Y_cluster.nrows || X_cluster.ncols != Y_cluster.ncols) eslog::error("size of X_cluster and Y_cluster does not match\n");
    if(X_cluster.order != Y_cluster.order) eslog::error("orders do not match\n");
    if(X_cluster.nrows != n_dofs_cluster_interface) eslog::error("incompatible cluster matrix size\n");
    char order = X_cluster.order;

    stacktimer::push("cpu_implicit_apply_matrix_compute");

    math::operations::fill_dnx<T>::do_all(&Y_cluster, T{0});

    _apply_hfeti_stuff(X_cluster_old, Y_cluster_old, filter);

    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & my = domain_data[di];

        VectorDenseView_new<I> my_D2C;
        my_D2C.set_view(feti.D2C[di].size(), feti.D2C[di].data(), AllocatorDummy_new::get_singleton(true, false));

        VectorDenseView_new<I> my_filter;
        my_filter.set_view(filter_transposed[di].size(), filter_transposed[di].data(), AllocatorDummy_new::get_singleton(true,false));
        size_t my_width = my_filter.size;

        MatrixDenseData_new<T> Z;
        Z.set(my.n_dofs_interface, my_width, order, AllocatorCPU_new::get_singleton());
        Z.alloc();

        MatrixDenseData_new<T> W;
        W.set(my.n_dofs_domain, my_width, order, AllocatorCPU_new::get_singleton());
        W.alloc();

        math::operations::submatrix_dnx_dnx_noncontig<T,I>::do_all(&X_cluster, &Z, &my_D2C, &my_filter);

        MatrixCsxView_new<T,I> B = MatrixCsxView_new<T,I>::from_old(feti.B1[di]);
        MatrixCsxView_new<T,I> Bt = B.get_transposed_reordered_view();

        math::operations::gemm_csx_dny_dnz<T,I>::do_all(&Bt, &Z, &W, T{1}, T{0});

        my.op_solver->solve(W, W);

        math::operations::gemm_csx_dny_dnz<T,I>::do_all(&B, &W, &Z, T{1}, T{0});

        math::operations::supermatrix_dnx_dnx_noncontig<T,I>::do_all(&Z, &Y_cluster, &my_D2C, &my_filter, math::operations::supermatrix_dnx_dnx_noncontig<T,I>::mode::accumulate_atomic);
    }
    if(!cfg.inner_timers) stacktimer::enable();

    stacktimer::pop();
    
    stacktimer::push("dual_synchronize");
    Y_cluster_old.synchronize();
    stacktimer::pop();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIImplicitGeneralSparseSolverCpu::toPrimal");
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
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::BtL(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("HybridFETIImplicitGeneralSparseSolverCpu::BtL");

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
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::_apply_hfeti_stuff(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    #pragma omp parallel for schedule(static,1)
    for (size_t di = 0; di < feti.K.size(); ++di) {
        applyBt(feti, di, x, Btx[di]);
    }

    _applyK(Btx, KplusBtx, false);

    applyB_threaded(feti, KplusBtx, y);
}



template <typename T, typename I>
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::_apply_hfeti_stuff(const Matrix_Dual<T> &x, Matrix_Dual<T> &y)
{
    Vector_Dual<T> _x, _y;
    _x.size = _y.size = x.ncols;

    for (int r = 0; r < x.nrows; ++r) {
        _x.vals = x.vals + x.ncols * r;
        _y.vals = y.vals + y.ncols * r;

        #pragma omp parallel for schedule(static,1)
        for (size_t di = 0; di < feti.K.size(); ++di) {
            applyBt(feti, di, _x, Btx[di]);
        }

        _applyK(Btx, KplusBtx, false);

        applyB_threaded(feti, KplusBtx, _y);
    }
}



template <typename T, typename I>
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::_apply_hfeti_stuff(const Matrix_Dual<T> &x, Matrix_Dual<T> &y, const std::vector<int> &filter)
{
    Vector_Dual<T> _x, _y;
    _x.size = _y.size = x.ncols;

    for (int r = 0; r < x.nrows; ++r) {
        if (filter[r]) {
            _x.vals = x.vals + x.ncols * r;
            _y.vals = y.vals + y.ncols * r;

            #pragma omp parallel for schedule(static,1)
            for (size_t di = 0; di < feti.K.size(); ++di) {
                applyBt(feti, di, _x, Btx[di]);
            }

            _applyK(Btx, KplusBtx, false);

            applyB_threaded(feti, KplusBtx, _y);
        }
    }
}



template <typename T, typename I>
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::_apply_hfeti_stuff(const Matrix_Dual<T> &x, Matrix_Dual<T> &y, const std::vector<std::vector<int>> &filter)
{
    eslog::error("not implemented\n");
}



// https://dl.acm.org/doi/pdf/10.1145/2929908.2929909

template <typename T, typename I>
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::_applyK(std::vector<Vector_Dense<T> > &b, std::vector<Vector_Dense<T> > &x, bool do_Kplus_solve)
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
            domain_data[di].op_solver->solve(rhs, sol);
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

        math::blas::multiply(T{-1}, dKB0[di], mu_di, T{1}, x[di], true);
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
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::_compute_beta_mu(std::vector<Vector_Dense<T> > &b)
{
    // g = B0 * (K+)^-1 * b
    math::set(g, T{0});
    #pragma omp parallel for schedule(static,1)
    for (size_t di = 0; di < feti.K.size(); ++di) {
        Vector_Dense<T> KB0b; KB0b.resize(dKB0[di].nrows);
        math::blas::apply(KB0b, T{1}, dKB0[di], T{0}, b[di]);
        for (size_t i = 0; i < D2C0[di].size(); ++i) {
            #pragma omp atomic
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
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::_computeB0()
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
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::_computeF0()
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
            domain_data[di].op_solver->solve(rhs, sol);
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
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::_computeG0()
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
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::_computeS0()
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
void HybridFETIImplicitGeneralSparseSolverCpu<T,I>::setup_config(config & cfg, const FETIConfiguration & feti_ecf_config)
{
    // defaults are set in config definition
    // if ecf value is auto, cfg value is not changed

    using ecf_config = DualopHybridfetiImplicitGeneralSparseSolverCpuConfig;
    const ecf_config & ecf = feti_ecf_config.dualop_hybridfeti_implicit_generalsparsesolver_cpu_config;

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
template class HybridFETIImplicitGeneralSparseSolverCpu<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

        #define INSTANTIATE \
        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        // INSTANTIATE_T(std::complex<double>)

            INSTANTIATE

    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
