
#include "totalfeti.explicit.sc.h"
#include "math/wrappers/math.blas.h"
#include "feti/common/applyB.h"
#include "basis/utilities/minmaxavg.h"
#include "my_timer.h"



namespace espreso {

template<typename T, typename I>
TotalFETIExplicitSc<T,I>::TotalFETIExplicitSc(FETI<T> &feti)
: DualOperator<T>(feti)
{
    if(!DirectSparseSolver<T,I>::provideSC()) {
        eslog::error("Cannot use explicit sc dual operator, sparse solver does not provide schur complement\n");
    }
}



template<typename T, typename I>
TotalFETIExplicitSc<T,I>::~TotalFETIExplicitSc()
{

}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::info()
{
    eslog::info(" = EXPLICIT TOTAL FETI OPERATOR USING SPBLAS SCHUR COMPLEMENT                                = \n");
    eslog::info(" =   EXTERNAL SCHUL COMPLEMENT SOLVER     %50s = \n", DirectSparseSolver<T,I>::name());
    eslog::info(minmaxavg<double>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.F.nrows * data.F.get_ld() * sizeof(T) / (1024.0 * 1024.0); }).to_string("  F MEMORY [MB]").c_str());
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_domain; }).to_string("  Domain volume [dofs]").c_str());
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_interface; }).to_string("  Domain surface [dofs]").c_str());
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::set(const step::Step &step)
{
    n_domains = feti.K.size();
    domain_data.resize(n_domains);
    for(size_t di = 0; di < n_domains; di++) {
        domain_data[di].n_dofs_interface = feti.B1[di].nrows;
        domain_data[di].n_dofs_domain = feti.B1[di].ncols;
    }

    constexpr bool Fs_share_memory = true;

    if constexpr(Fs_share_memory) {
        Fs_allocated.resize((n_domains - 1) / 2 + 1);
        std::vector<size_t> domain_idxs_sorted_by_f_size_desc(n_domains);
        for(size_t di = 0; di < n_domains; di++) {
            domain_idxs_sorted_by_f_size_desc[di] = di;
        }
        std::sort(domain_idxs_sorted_by_f_size_desc.rbegin(), domain_idxs_sorted_by_f_size_desc.rend(), [&](size_t dl, size_t dr){ return domain_data[dl].n_dofs_interface < domain_data[dr].n_dofs_interface; });
        for(size_t i = 0; i < n_domains; i++) {
            size_t di = domain_idxs_sorted_by_f_size_desc[i];
            auto & data = domain_data[di];

            if(i % 2 == 0) {
                Fs_allocated[i / 2].resize(data.n_dofs_interface + 1, data.n_dofs_interface);
                data.F_fill = 'U';
                data.F.shallowCopy(Fs_allocated[i / 2]);
                data.F.nrows = data.n_dofs_interface;
                data.F.ncols = data.n_dofs_interface;
            }
            else {
                data.F_fill = 'L';
                data.F.shallowCopy(Fs_allocated[i / 2]);
                data.F.nrows = data.n_dofs_interface;
                data.F.ncols = data.n_dofs_interface;
                data.F.vals += data.F.get_ld();
            }
            if constexpr(utils::is_real<T>())    data.F.type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE;
            if constexpr(utils::is_complex<T>()) data.F.type = Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE;
        }
    }
    else {
        for(size_t di = 0; di < n_domains; di++) {
            auto & data = domain_data[di];

            data.F_fill = 'U';
            data.F.resize(data.n_dofs_interface, data.n_dofs_interface, data.n_dofs_interface);
            if constexpr(utils::is_real<T>())    data.F.type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE;
            if constexpr(utils::is_complex<T>()) data.F.type = Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE;
        }
    }

    my_timer tm_total, tm_inner, tm_trans, tm_commit, tm_symbfact;

    tm_total.start();
    #pragma omp parallel for
    for(size_t di = 0; di < n_domains; di++) {
        tm_inner.start();
        auto & data = domain_data[di];

        if(getSymmetry(feti.K[di].type) != Matrix_Symmetry::HERMITIAN || feti.K[di].shape != Matrix_Shape::UPPER) {
            eslog::error("implemented only for hermitian K matrices stored in upper triangle. TODO\n");
        }

        data.Kreg.type = feti.K[di].type;
        data.Kreg.shape = feti.K[di].shape;
        math::combine(data.Kreg, feti.K[di], feti.RegMat[di]);

        Matrix_CSR<T,I> & B = feti.B1[di];
        data.Bt.resize(B.ncols, B.nrows, B.nnz);
        if constexpr(utils::is_real<T>())    data.Bt.type = Matrix_Type::REAL_NONSYMMETRIC;
        if constexpr(utils::is_complex<T>()) data.Bt.type = Matrix_Type::COMPLEX_NONSYMMETRIC;
        data.map_B_transpose.resize(B.nnz);
        math::csrTranspose(data.Bt, B, data.map_B_transpose, math::CsrTransposeStage::Pattern, true);

        data.null_matrix_A21.resize(data.n_dofs_interface, data.n_dofs_domain, 0);
        std::fill_n(data.null_matrix_A21.rows, data.null_matrix_A21.nrows + 1, 0);
        if constexpr(utils::is_real<T>())    data.null_matrix_A21.type = Matrix_Type::REAL_NONSYMMETRIC;
        if constexpr(utils::is_complex<T>()) data.null_matrix_A21.type = Matrix_Type::COMPLEX_NONSYMMETRIC;
        data.null_matrix_A22.resize(data.n_dofs_interface, data.n_dofs_interface, 0);
        std::fill_n(data.null_matrix_A22.rows, data.null_matrix_A22.nrows + 1, 0);
        if constexpr(utils::is_real<T>())    data.null_matrix_A22.type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE;
        if constexpr(utils::is_complex<T>()) data.null_matrix_A22.type = Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE;
        data.null_matrix_A22.shape = Matrix_Shape::UPPER;

        tm_commit.start();
        data.sc_solver.commitMatrix(data.Kreg, data.Bt, data.null_matrix_A21, data.null_matrix_A22);
        tm_commit.stop();
        tm_symbfact.start();
        data.sc_solver.factorizeSymbolic();
        tm_symbfact.stop();

        data.x.resize(data.n_dofs_interface);
        data.y.resize(data.n_dofs_interface);
        tm_inner.stop();
    }
    tm_total.stop();

    print_timer("Set     total", tm_total);
    print_timer("Set       inner", tm_inner);
    print_timer("Set         commit", tm_commit);
    print_timer("Set         symbfact", tm_symbfact);
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::update(const step::Step &step)
{
    my_timer tm_total, tm_inner, tm_updatevals, tm_factnumer, tm_maked;

    tm_total.start();
    #pragma omp parallel for
    for(size_t di = 0; di < n_domains; di++) {
        tm_inner.start();
        auto & data = domain_data[di];

        math::sumCombined(data.Kreg, T{1.0}, feti.K[di], feti.RegMat[di]);

        Matrix_CSR<T,I> & B = feti.B1[di];
        math::csrTranspose(data.Bt, B, data.map_B_transpose, math::CsrTransposeStage::Values, true);

        tm_updatevals.start();
        data.sc_solver.updateMatrixValues();
        tm_updatevals.stop();
        tm_factnumer.start();
        data.sc_solver.factorizeNumericAndGetSc(data.F, data.F_fill, T{-1});
        tm_factnumer.stop();
        tm_inner.stop();
    }
    tm_total.stop();

    tm_maked.start();
    {
        if (feti.updated.B) {
            d.resize();
        }
        std::vector<Vector_Dense<T,I>> Kplus_fs(n_domains);
        #pragma omp parallel for schedule(static,1)
        for(size_t di = 0; di < n_domains; di++) {
            auto & data = domain_data[di];
            Kplus_fs[di].resize(data.n_dofs_domain);
            data.sc_solver.solveA11(feti.f[di], Kplus_fs[di]);
        }
        applyB(feti, Kplus_fs, d);
        math::add(d, T{-1}, feti.c);
    }
    tm_maked.stop();

    print_timer("Update  total", tm_total);
    print_timer("Update    inner", tm_inner);
    print_timer("Update      updatevals", tm_updatevals);
    print_timer("Update      factnumer", tm_factnumer);
    print_timer("Update  maked", tm_maked);
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::_apply(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    my_timer tm_total, tm_mv;

    double start = eslog::time();

    tm_total.start();
    memset(y_cluster.vals, 0, y_cluster.size * sizeof(T));

    #pragma omp parallel for
    for(size_t di = 0; di < n_domains; di++) {
        auto & data = domain_data[di];

        std::vector<I> & D2C = feti.D2C[di];

        for(I i = 0; i < data.n_dofs_interface; i++) {
            data.x.vals[i] = x_cluster.vals[D2C[i]];
        }

        tm_mv.start();
        math::blas::apply_hermitian(data.y, T{1}, data.F, data.F_fill, T{0}, data.x);
        tm_mv.stop();

        for(I i = 0; i < data.n_dofs_interface; i++) {
            #pragma omp atomic
            y_cluster.vals[D2C[i]] += data.y.vals[i];
        }
    }
    tm_total.stop();

    double stop = eslog::time();
    eslog::info("TMP DUAL OPERATOR APPLY TIME:  %12.6f ms\n", (stop - start) * 1000.0);

    y_cluster.synchronize();

    print_timer("Apply   total", tm_total);
    print_timer("Apply     mv", tm_mv);
}



template <typename T, typename I>
void TotalFETIExplicitSc<T,I>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    _apply(x, y);
    y.synchronize();
}



template <typename T, typename I>
void TotalFETIExplicitSc<T,I>::apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y)
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
void TotalFETIExplicitSc<T,I>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    #pragma omp parallel for schedule(static,1)
    for (size_t di = 0; di < n_domains; ++di) {
        Vector_Dense<T,I> z;
        z.resize(y[di]);
        applyBt(feti, di, x, z, T{-1});
        math::add(z, T{1}, feti.f[di]);
        domain_data[di].sc_solver.solveA11(z, y[di]);
    }
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::print(const step::Step &step)
{
    eslog::error("TotalFETIExplicitSc::print not implemented");
}





#define INSTANTIATE_T_I(T,I) \
template class TotalFETIExplicitSc<T,I>;

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
