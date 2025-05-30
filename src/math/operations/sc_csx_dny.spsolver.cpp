
#include "math/operations/sc_csx_dny.spsolver.h"

#include "math/primitives_new/allocator_new.h"
#include "math/wrappers/math.spsolver.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/matrix_csx_data_new.h"
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
struct sc_csx_dny_spsolver_data
{
    DirectSparseSolver<T,I> A11_solver;
    Matrix_CSR<T,I> A11_old;
    quadrisect_csx_csy<T,I> op_split;
    MatrixCsxData_new<T,I> sub_A11;
    MatrixCsxData_new<T,I> sub_A12;
    MatrixCsxData_new<T,I> sub_A21;
    MatrixCsxData_new<T,I> sub_A22;
};



template<typename T, typename I>
sc_csx_dny_spsolver<T,I>::sc_csx_dny_spsolver() = default;



template<typename T, typename I>
sc_csx_dny_spsolver<T,I>::~sc_csx_dny_spsolver() = default;



template<typename T, typename I>
void sc_csx_dny_spsolver<T,I>::internal_preprocess()
{
    data = std::make_unique<sc_csx_dny_spsolver_data<T,I>>();

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

    data->A11_old = MatrixCsxView_new<T,I>::template to_old<cpu_allocator>(*A11);

    data->A11_solver.commit(data->A11_old);
    data->A11_solver.symbolicFactorization();
}



template<typename T, typename I>
void sc_csx_dny_spsolver<T,I>::internal_perform_1()
{
    if(called_set_matrix == '1') {
        data->op_split.perform();
    }

    data->A11_old = MatrixCsxView_new<T,I>::template to_old<cpu_allocator>(*A11);

    data->A11_solver.numericalFactorization();
}



template<typename T, typename I>
void sc_csx_dny_spsolver<T,I>::internal_perform_2()
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

    Matrix_Dense<T,I> X_old(size_sc, size_A11);
    MatrixDenseView_new<T> X = MatrixDenseView_new<T>::from_old(X_old, 'R').get_transposed_reordered_view();
    X.prop.uplo = 'F';

    Matrix_Dense<T,I> Y_old(size_sc, size_A11);
    MatrixDenseView_new<T> Y = MatrixDenseView_new<T>::from_old(Y_old, 'R').get_transposed_reordered_view();
    Y.prop.uplo = 'F';

    math::operations::convert_csx_dny<T,I>::do_all(A12_to_use, &X);

    data->A11_solver.solve(X_old, Y_old);

    if(is_matrix_hermitian) {
        // carefull not to touch the other triangle of sc

        MatrixDenseData_new<T> sc_tmp;
        sc_tmp.set(size_sc, size_sc, sc->order, AllocatorCPU_new::get_singleton());
        sc_tmp.alloc();
        sc_tmp.prop.uplo = sc->prop.uplo;

        if(A22 == nullptr) math::operations::fill_dnx<T>::do_all(&sc_tmp, T{0});
        if(A22 != nullptr) math::operations::convert_csx_dny<T,I>::do_all(A22, &sc_tmp);

        math::operations::gemm_csx_dny_dnz<T,I>::do_all(A21_to_use, &Y, &sc_tmp, T{-1} * alpha, T{1} * alpha);

        math::operations::copy_dnx<T>::do_all(&sc_tmp, sc);
    }
    else {
        if(A22 == nullptr) math::operations::fill_dnx<T>::do_all(sc, T{0});
        if(A22 != nullptr) math::operations::convert_csx_dny<T,I>::do_all(A22, sc);

        math::operations::gemm_csx_dny_dnz<T,I>::do_all(A21_to_use, &Y, sc, T{-1} * alpha, T{1} * alpha);
    }
}



template<typename T, typename I>
void sc_csx_dny_spsolver<T,I>::internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    Vector_Dense<T,I> rhs_old = VectorDenseView_new<T>::template to_old<I,cpu_allocator>(rhs);
    Vector_Dense<T,I> sol_old = VectorDenseView_new<T>::template to_old<I,cpu_allocator>(sol);

    data->A11_solver.solve(rhs_old, sol_old);
}



#define INSTANTIATE_T_I(T,I) \
template class sc_csx_dny_spsolver<T,I>;

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
