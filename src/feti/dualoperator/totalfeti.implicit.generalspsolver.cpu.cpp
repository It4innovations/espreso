
#include "feti/dualoperator/totalfeti.implicit.generalspsolver.cpu.h"

#include "math/primitives_new/allocator_new.h"
#include "feti/common/applyB.h"
#include "math/math.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/minmaxavg.h"
#include "basis/utilities/stacktimer.h"
#include "math/operations/gemv_csx.h"



namespace espreso {



template<typename T, typename I>
TotalFETIImplicitGeneralSpSolverCpu<T,I>::TotalFETIImplicitGeneralSpSolverCpu(FETI<T> &feti) : DualOperator<T>(feti)
{
    
}



template<typename T, typename I>
TotalFETIImplicitGeneralSpSolverCpu<T,I>::~TotalFETIImplicitGeneralSpSolverCpu()
{
}



template<typename T, typename I>
void TotalFETIImplicitGeneralSpSolverCpu<T,I>::info()
{
    eslog::info(" = IMPLICIT TOTAL FETI OPERATOR USING GENERAL SPARSE SOLVER ON CPU                    = \n");

    if(cfg.print_parameters) {
        auto bool_to_string = [](bool val){ return val ? "TRUE" : "FALSE";};

        eslog::info(" =   %-50s       %+30s = \n", "parallel_set", bool_to_string(cfg.parallel_set));
        eslog::info(" =   %-50s       %+30s = \n", "parallel_update", bool_to_string(cfg.parallel_update));
        eslog::info(" =   %-50s       %+30s = \n", "parallel_apply", bool_to_string(cfg.parallel_apply));
    }

    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_domain; }).to_string("  Domain volume [dofs]").c_str());
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_interface; }).to_string("  Domain surface [dofs]").c_str());
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}



template<typename T, typename I>
void TotalFETIImplicitGeneralSpSolverCpu<T,I>::set(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIImplicitGeneralSpSolverCpu::set");

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

        my.op_solver = math::operations::solver_csx<T,I>::make(cfg.solver_is, &my.Kreg, false, true);
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
void TotalFETIImplicitGeneralSpSolverCpu<T,I>::update(const step::Step &step)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIImplicitGeneralSpSolverCpu::update");

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
void TotalFETIImplicitGeneralSpSolverCpu<T,I>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIImplicitGeneralSpSolverCpu::apply");

    std::fill_n(y.vals, y.size, T{0});

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

    y.synchronize();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void TotalFETIImplicitGeneralSpSolverCpu<T,I>::apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y)
{
    eslog::error("not implemented yet\n");
    // Vector_Dual<T> _x, _y;
    // _x.size = _y.size = x.ncols;
    // for (int r = 0; r < x.nrows; ++r) {
    //     _x.vals = x.vals + x.ncols * r;
    //     _y.vals = y.vals + y.ncols * r;
    //     _apply(_x, _y);
    // }
    // y.synchronize();
}



template<typename T, typename I>
void TotalFETIImplicitGeneralSpSolverCpu<T,I>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("TotalFETIImplicitGeneralSpSolverCpu::toPrimal");
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
void TotalFETIImplicitGeneralSpSolverCpu<T,I>::print(const step::Step &step)
{
    eslog::error("TotalFETIImplicitGeneralSpSolverCpu::print not implemented");
}



#define INSTANTIATE_T_I(T,I) \
template class TotalFETIImplicitGeneralSpSolverCpu<T,I>;

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
