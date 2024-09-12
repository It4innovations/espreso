
#include "totalfeti.explicit.sc.h"
#include "math/wrappers/math.blas.h"
#include "feti/common/applyB.h"



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
    // TODO
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::set(const step::Step &step)
{
    n_domains = feti.K.size();
    domain_data.resize(n_domains);

    // todo share memory of F triangles

    #pragma omp parallel for schedule(static,1)
    for(size_t di = 0; di < n_domains; di++) {
        auto & data = domain_data[di];

        if(getSymmetry(feti.K[di].type) != Matrix_Symmetry::HERMITIAN || feti.K[di].shape != Matrix_Shape::UPPER) {
            eslog::error("implemented only for hermitian K matrices stored in upper triangle\n");
        }

        data.n_dofs_interface = feti.B1[di].nrows;
        data.n_dofs_domain = feti.B1[di].ncols;

        data.Kreg.type = feti.K[di].type;
        data.Kreg.shape = feti.K[di].shape;
        math::combine(data.Kreg, feti.K[di], feti.RegMat[di]);

        Matrix_CSR<T,I> & B = feti.B1[di];
        data.Bt.resize(B.ncols, B.nrows, B.nnz);
        data.map_B_transpose.resize(B.nnz);
        math::csrTranspose(data.Bt, B, data.map_B_transpose, math::CsrTransposeStage::Pattern, true);

        data.null_matrix_A21.resize(data.n_dofs_interface, data.n_dofs_domain, 0);
        std::fill_n(data.null_matrix_A21.rows, data.null_matrix_A21.nrows + 1, 0);
        data.null_matrix_A22.resize(data.n_dofs_interface, data.n_dofs_interface, 0);
        std::fill_n(data.null_matrix_A22.rows, data.null_matrix_A22.nrows + 1, 0);

        data.sc_solver.commitMatrix(data.Kreg, data.Bt, data.null_matrix_A21, data.null_matrix_A22);
        data.sc_solver.factorizeSymbolic();

        data.F.resize(data.n_dofs_interface, data.n_dofs_interface);
        if constexpr(utils::is_real<T>())    data.F.type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE;
        if constexpr(utils::is_complex<T>()) data.F.type = Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE;

        data.x.resize(data.n_dofs_interface);
        data.y.resize(data.n_dofs_interface);
    }
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::update(const step::Step &step)
{
    #pragma omp parallel for schedule(static,1)
    for(size_t di = 0; di < n_domains; di++) {
        auto & data = domain_data[di];

        // when new SCsolver will exist, have something like prefered_fill (uplo), preferred_input (one big or 4 small matrices), etc.

        math::sumCombined(data.Kreg, T{1.0}, feti.K[di], feti.RegMat[di]);

        Matrix_CSR<T,I> & B = feti.B1[di];
        math::csrTranspose(data.Bt, B, data.map_B_transpose, math::CsrTransposeStage::Values, true);

        data.sc_solver.updateMatrixValues();
        data.sc_solver.factorizeNumericAndGetSc(data.F, 'U');

        math::blas::scale(data.F.nrows * data.F.get_ld(), T{-1}, data.F.vals, 1);
    }
    
    {
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
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::_apply(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    memset(y_cluster.vals, 0, y_cluster.size * sizeof(T));

    #pragma omp parallel for schedule(static,1)
    for(size_t di = 0; di < n_domains; di++) {
        auto & data = domain_data[di];

        std::vector<I> & D2C = feti.D2C[di];

        for(I i = 0; i < data.n_dofs_interface; i++) {
            data.x.vals[i] = x_cluster.vals[D2C[i]];
        }

        math::blas::apply_hermitian(data.y, T{1}, data.F, 'U', T{0}, data.x);

        for(I i = 0; i < data.n_dofs_interface; i++) {
            #pragma omp atomic
            y_cluster.vals[D2C[i]] += data.y.vals[i];
        }
    }

    y_cluster.synchronize();
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
