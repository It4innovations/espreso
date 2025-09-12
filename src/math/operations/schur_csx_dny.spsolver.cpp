
#include "math/operations/schur_csx_dny.spsolver.h"

#include "config/ecf/operations/schur_csx_dny.spsolver.h"
#include "esinfo/ecfinfo.h"
#include "math/primitives_new/allocator_new.h"
#include "math/operations/solver_csx.h"
#include "math/operations/convert_csx_dny.h"
#include "math/operations/convert_dnx_dny.h"
#include "math/operations/copy_dnx.h"
#include "math/operations/fill_dnx.h"
#include "math/operations/quadrisect_csx_csy.h"
#include "math/operations/gemm_csx_dny_dnz.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
struct schur_csx_dny_spsolver_data
{
    struct config
    {
        typename solver_csx<T,I>::implementation_selector a11_solver_impl = solver_csx<T,I>::implementation_selector::autoselect;
        char X_order = 'C';
    } cfg;
    std::unique_ptr<solver_csx<T,I>> op_A11_solver;
    quadrisect_csx_csy<T,I> op_split;
    MatrixCsxData_new<T,I> sub_A11;
    MatrixCsxData_new<T,I> sub_A12;
    MatrixCsxData_new<T,I> sub_A21;
    MatrixCsxData_new<T,I> sub_A22;
};



template<typename T, typename I>
void setup_config(typename schur_csx_dny_spsolver_data<T,I>::config & cfg)
{
    using solver_impl_t = typename solver_csx<T,I>::implementation_selector;
    using ecf_config = SchurCsxDnySpsolverConfig;
    const ecf_config & ecf = info::ecf->operations.schur_csx_dny_spsolver;

    switch(ecf.order_X) {
        case ecf_config::MATRIX_ORDER::AUTO: break;
        case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.X_order = 'R'; break;
        case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.X_order = 'C'; break;
    }

    switch(ecf.sparse_solver_impl) {
        case ecf_config::SPARSE_SOLVER_IMPL::AUTO: break;
        case ecf_config::SPARSE_SOLVER_IMPL::MKLPARDISO:   cfg.a11_solver_impl = solver_impl_t::mklpardiso;   break;
        case ecf_config::SPARSE_SOLVER_IMPL::SUITESPARSE:  cfg.a11_solver_impl = solver_impl_t::suitesparse;  break;
        case ecf_config::SPARSE_SOLVER_IMPL::MUMPS:        cfg.a11_solver_impl = solver_impl_t::mumps;        break;
        case ecf_config::SPARSE_SOLVER_IMPL::STRUMPACK:    cfg.a11_solver_impl = solver_impl_t::strumpack;    break;
        case ecf_config::SPARSE_SOLVER_IMPL::PASTIX:       cfg.a11_solver_impl = solver_impl_t::pastix;       break;
        case ecf_config::SPARSE_SOLVER_IMPL::SUPERLU_DIST: cfg.a11_solver_impl = solver_impl_t::superlu_dist; break;
    }
}



template<typename T, typename I>
schur_csx_dny_spsolver<T,I>::schur_csx_dny_spsolver() {}



template<typename T, typename I>
schur_csx_dny_spsolver<T,I>::~schur_csx_dny_spsolver() {}



template<typename T, typename I>
void schur_csx_dny_spsolver<T,I>::internal_preprocess()
{
    data = std::make_unique<schur_csx_dny_spsolver_data<T,I>>();
    setup_config<T,I>(data->cfg);

    if(called_set_matrix == '1') {
        data->op_split.set_matrix_src(A);
        data->op_split.set_matrices_dst(&data->sub_A11, &data->sub_A12, &data->sub_A21, &data->sub_A22);
        data->op_split.set_bounds(size_A11, size_A11);
        data->op_split.setup();

        data->sub_A11.set(size_A11, size_A11, data->op_split.get_output_matrix_11_nnz(), A->order, AllocatorCPU_new::get_singleton());
        data->sub_A12.set(size_A11, size_sc,  data->op_split.get_output_matrix_12_nnz(), A->order, AllocatorCPU_new::get_singleton());
        data->sub_A21.set(size_sc,  size_A11, data->op_split.get_output_matrix_21_nnz(), A->order, AllocatorCPU_new::get_singleton());
        data->sub_A22.set(size_sc,  size_sc,  data->op_split.get_output_matrix_22_nnz(), A->order, AllocatorCPU_new::get_singleton());

        data->sub_A11.prop = A->prop;
        data->sub_A22.prop = A->prop;

        data->sub_A11.alloc();
        data->sub_A12.alloc();
        data->sub_A21.alloc();
        data->sub_A22.alloc();

        data->op_split.perform();

        A11 = &data->sub_A11;
        if(A->prop.uplo != 'L') A12 = &data->sub_A12;
        if(A->prop.uplo != 'U') A21 = &data->sub_A21;
        A22 = &data->sub_A22;
    }

    data->op_A11_solver = solver_csx<T,I>::make(data->cfg.a11_solver_impl, &A11->prop, false, true);
    data->op_A11_solver->set_matrix_A(A11);
    data->op_A11_solver->set_needs(false, true);
    data->op_A11_solver->factorize_symbolic();
}



template<typename T, typename I>
void schur_csx_dny_spsolver<T,I>::internal_perform_1()
{
    if(called_set_matrix == '1') {
        data->op_split.perform();
    }

    data->op_A11_solver->factorize_numeric();
}



template<typename T, typename I>
void schur_csx_dny_spsolver<T,I>::internal_perform_2()
{
    MatrixCsxView_new<T,I> * A12_to_use = A12;
    MatrixCsxView_new<T,I> A21_transposed_reordered;
    if(A12 == nullptr && is_matrix_hermitian) {
        A21_transposed_reordered = A21->get_transposed_reordered_view();
        A12_to_use = &A21_transposed_reordered;
    }

    MatrixCsxView_new<T,I> * A21_to_use = A21;
    MatrixCsxView_new<T,I> A12_transposed_reordered;
    if(A21 == nullptr && is_matrix_hermitian) {
        A12_transposed_reordered = A12->get_transposed_reordered_view();
        A21_to_use = &A12_transposed_reordered;
    }

    MatrixDenseData_new<T> X;
    X.set(size_A11, size_sc, data->cfg.X_order, AllocatorCPU_new::get_singleton());
    X.alloc();

    data->op_A11_solver->solve(*A12_to_use, X);

    if(is_matrix_hermitian) {
        // carefull not to touch the other triangle of sc

        MatrixDenseData_new<T> sc_tmp;
        sc_tmp.set(size_sc, size_sc, sc->order, AllocatorCPU_new::get_singleton());
        sc_tmp.alloc();
        sc_tmp.prop.uplo = sc->prop.uplo;

        if(A22 == nullptr) {
            math::operations::fill_dnx<T>::do_all(&sc_tmp, T{0});
        }
        if(A22 != nullptr) {
            MatrixCsxView_new<T,I> A22_rt = A22->get_transposed_reordered_view();
            MatrixCsxView_new<T,I> * A22_to_use = (A22->prop.uplo == sc_tmp.prop.uplo) ? A22 : &A22_rt;
            math::operations::convert_csx_dny<T,I>::do_all(A22_to_use, &sc_tmp);
        }

        math::operations::gemm_csx_dny_dnz<T,I>::do_all(A21_to_use, &X, &sc_tmp, T{-1} * alpha, T{1} * alpha);

        math::operations::copy_dnx<T>::do_all(&sc_tmp, sc);
    }
    else {
        if(A22 == nullptr) math::operations::fill_dnx<T>::do_all(sc, T{0});
        if(A22 != nullptr) math::operations::convert_csx_dny<T,I>::do_all(A22, sc);

        math::operations::gemm_csx_dny_dnz<T,I>::do_all(A21_to_use, &X, sc, T{-1} * alpha, T{1} * alpha);
    }
}



template<typename T, typename I>
void schur_csx_dny_spsolver<T,I>::internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    data->op_A11_solver->solve(rhs, sol);
}



#define INSTANTIATE_T_I(T,I) \
template class schur_csx_dny_spsolver<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T,int32_t) \
    /* INSTANTIATE_T_I(T,int64_t) */

        #define INSTANTIATE \
        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        INSTANTIATE_T(std::complex<double>)

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}
