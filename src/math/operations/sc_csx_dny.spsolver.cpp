
#include "math/operations/sc_csx_dny.spsolver.h"

#include "math/primitives_new/allocator_new.h"
#include "math/wrappers/math.spsolver.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/operations/convert_csx_dny.h"
#include "math/operations/convert_dnx_dny.h"
#include "math/operations/copy_dnx.h"
#include "math/operations/fill_dnx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
struct sc_csx_dny_spsolver_data
{
    DirectSparseSolver<T,I> A11_solver;
    Matrix_CSR<T,I> A11_old;
};



template<typename T, typename I>
sc_csx_dny_spsolver<T,I>::sc_csx_dny_spsolver() = default;



template<typename T, typename I>
sc_csx_dny_spsolver<T,I>::~sc_csx_dny_spsolver() = default;



template<typename T, typename I>
void sc_csx_dny_spsolver<T,I>::internal_preprocess()
{
    if(called_set_matrix == '1') eslog::error("dont yet support single large input matrix\n");

    data = std::make_unique<sc_csx_dny_spsolver_data<T,I>>();

    data->A11_old = MatrixCsxView_new<T,I>::template to_old<cpu_allocator>(*A11);

    data->A11_solver.commit(data->A11_old);
    data->A11_solver.symbolicFactorization();
}



template<typename T, typename I>
void sc_csx_dny_spsolver<T,I>::internal_perform_1()
{
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

    // MatrixDenseData_new<T> X;
    // X.set(size_A11, size_sc, 'C', AllocatorCPU_new::get_singleton());
    // X.alloc();
    // X.prop.uplo = 'F';
    // Matrix_Dense<T,I> X_old = MatrixDenseView_new<T>::template to_old<I,cpu_allocator>(X);
    
    // MatrixDenseData_new<T> Y;
    // Y.set(size_A11, size_sc, 'C', AllocatorCPU_new::get_singleton());
    // Y.alloc();
    // Y.prop.uplo = 'F';
    // Matrix_Dense<T,I> Y_old = MatrixDenseView_new<T>::template to_old<I,cpu_allocator>(Y);

    Matrix_Dense<T,I> X_old(size_sc, size_A11);
    MatrixDenseView_new<T> X = MatrixDenseView_new<T>::from_old(X_old, 'R').get_transposed_reordered_view();
    X.prop.uplo = 'F';

    Matrix_Dense<T,I> Y_old(size_sc, size_A11);
    MatrixDenseView_new<T> Y = MatrixDenseView_new<T>::from_old(Y_old, 'R').get_transposed_reordered_view();
    Y.prop.uplo = 'F';

    math::operations::convert_csx_dny<T,I>::do_all(A12_to_use, &X);

    data->A11_solver.solve(X_old, Y_old);

    MatrixDenseData_new<T> sc_tmp;
    sc_tmp.set(size_sc, size_sc, sc->order, AllocatorCPU_new::get_singleton());
    sc_tmp.alloc();
    sc_tmp.prop.uplo = sc->prop.uplo;

    if(Y.order == sc_tmp.order) {
        if(A22 == nullptr) math::operations::fill_dnx<T>::do_all(&sc_tmp, T{0});
        if(A22 != nullptr) math::operations::convert_csx_dny<T,I>::do_all(A22, &sc_tmp);

        math::spblas::handle_mm handlemm;
        math::spblas::mm<T,I>(*A21_to_use, Y, sc_tmp, T{-1}, T{1}, handlemm, 'A');
    }
    else {
        MatrixDenseData_new<T> sc_tmp2;
        sc_tmp2.set(size_sc, size_sc, Y.order, AllocatorCPU_new::get_singleton());
        sc_tmp2.alloc();
        sc_tmp2.prop.uplo = sc->prop.uplo;
        
        if(A22 == nullptr) math::operations::fill_dnx<T>::do_all(&sc_tmp2, T{0});
        if(A22 != nullptr) math::operations::convert_csx_dny<T,I>::do_all(A22, &sc_tmp2);

        math::spblas::handle_mm handlemm;
        math::spblas::mm<T,I>(*A21_to_use, Y, sc_tmp2, T{-1}, T{1}, handlemm, 'A');

        convert_dnx_dny<T>::do_all(&sc_tmp2, &sc_tmp, false);
    }

    math::operations::copy_dnx<T>::do_all(&sc_tmp, sc);
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
