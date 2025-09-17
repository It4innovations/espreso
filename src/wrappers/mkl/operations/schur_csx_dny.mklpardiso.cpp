
#ifdef HAVE_MKL

#include "wrappers/mkl/operations/schur_csx_dny.mklpardiso.h"

#include <mkl.h>

#include "math/primitives_new/allocator_new.h"
#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/vector_dense_data_new.h"
#include "math/operations/concat_csx.h"
#include "math/operations/lincomb_matrix_dnx.h"
#include "math/operations/convert_dnx_dny.h"
#include "math/operations/fill_dnx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
struct schur_csx_dny_mklpardiso_data
{
    MatrixCsxData_new<T,I> A_whole;
    VectorDenseData_new<MKL_INT> perm;
    VectorDenseData_new<T> x;
    VectorDenseData_new<T> y;
    void * pt[64];
    MKL_INT iparm[64];
    MKL_INT mtype;
    MKL_INT msglvl = 0;
    MKL_INT error;
    MKL_INT size_matrix_mklint;
};



template<typename T, typename I>
static int get_pardisoType(const MatrixCsxView_new<T,I> & M)
{
    Matrix_Type mt_old = get_old_matrix_type<T>(M);
    switch (mt_old) {
    case Matrix_Type::UNSET_INVALID_NONE:                    eslog::error("Invalid/unset matrix type\n");
    case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:    return  2;
    case Matrix_Type::REAL_SYMMETRIC_INDEFINITE:           return -2;
    case Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC:         return  1;
    case Matrix_Type::REAL_NONSYMMETRIC:                   return 11;
    case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE: return  4;
    case Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE:        return -4;
    case Matrix_Type::COMPLEX_SYMMETRIC:                   return  6;
    case Matrix_Type::COMPLEX_STRUCTURALLY_SYMMETRIC:      return  3;
    case Matrix_Type::COMPLEX_NONSYMMETRIC:                return 13;
    }
    return 0;
}



template<typename T, typename I>
schur_csx_dny_mklpardiso<T,I>::schur_csx_dny_mklpardiso() {}



template<typename T, typename I>
schur_csx_dny_mklpardiso<T,I>::~schur_csx_dny_mklpardiso()
{
    MKL_INT phase = -1;
    pardiso(data->pt, nullptr, nullptr, &data->mtype, &phase, &data->size_matrix_mklint, nullptr, nullptr, nullptr, nullptr, nullptr, data->iparm, &data->msglvl, nullptr, nullptr, &data->error);
    if(data->error != 0) {
        eslog::error("pardiso error %d\n", (int)data->error);
    }
}



template<typename T, typename I>
void schur_csx_dny_mklpardiso<T,I>::internal_preprocess()
{
    if(called_set_matrix == '1' && A->order != 'R') eslog::error("single input matrix needs to be rowmajor\n");

    data = std::make_unique<schur_csx_dny_mklpardiso_data<T,I>>();

    if(size_matrix > (size_t)std::numeric_limits<MKL_INT>::max()) eslog::error("size_matrix too large for MKL_INT\n");
    data->size_matrix_mklint = (MKL_INT)size_matrix;

    if(called_set_matrix == '4') {
        this->helper_concat(data->A_whole, 'P');
    }

    std::fill_n(data->pt, 64, nullptr);
    std::fill_n(data->iparm, 64, 0);
    data->iparm[0] = 1; // I did popullate iparm
    data->iparm[1] = 2; // metis fill reduction
    data->iparm[7] = 0; // Max number of refinement iterations (default changed between 2024.2 and 2025.0)
    data->iparm[34] = 1; // zero-based indexing
    data->iparm[35] = 1; // schur complement

    data->mtype = get_pardisoType(*A);

    data->perm.set(size_matrix, AllocatorCPU_new::get_singleton());
    data->perm.alloc();
    std::fill_n(data->perm.vals, size_A11, MKL_INT{0});
    std::fill_n(data->perm.vals + size_A11, size_sc, MKL_INT{1});

    MKL_INT phase = 11;
    MKL_INT one = 1;
    pardiso(data->pt, &one, &one, &data->mtype, &phase, &data->size_matrix_mklint, A->vals, A->ptrs, A->idxs, data->perm.vals, &one, data->iparm, &data->msglvl, nullptr, nullptr, &data->error);
    if(data->error != 0) {
        eslog::error("pardiso error %d\n", (int)data->error);
    }

    if(need_solve_A11) {
        data->x.set(size_matrix, AllocatorCPU_new::get_singleton());
        data->y.set(size_matrix, AllocatorCPU_new::get_singleton());
        data->x.alloc();
        data->y.alloc();
    }
}



template<typename T, typename I>
void schur_csx_dny_mklpardiso<T,I>::internal_perform_1()
{
    if(called_set_matrix == '4') {
        this->helper_concat(data->A_whole, 'F');
    }

    T * sc_tmp_vals = new T[size_sc * size_sc];
    MatrixDenseView_new<T> sc_tmp;
    sc_tmp.set_view(size_sc, size_sc, size_sc, sc->order, sc_tmp_vals, AllocatorCPU_new::get_singleton());

    MKL_INT phase = 22;
    MKL_INT one = 1;
    pardiso(data->pt, &one, &one, &data->mtype, &phase, &data->size_matrix_mklint, A->vals, A->ptrs, A->idxs, data->perm.vals, &one, data->iparm, &data->msglvl, nullptr, sc_tmp.vals, &data->error);
    if(data->error != 0) {
        eslog::error("pardiso error %d\n", (int)data->error);
    }

    lincomb_matrix_dnx<T>::do_all(sc, alpha, &sc_tmp, 0, nullptr);

    delete[] sc_tmp_vals;
}



template<typename T, typename I>
void schur_csx_dny_mklpardiso<T,I>::internal_perform_2()
{
}



template<typename T, typename I>
void schur_csx_dny_mklpardiso<T,I>::internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    MKL_INT phase;
    MKL_INT one = 1;

    std::copy_n(rhs.vals, size_A11, data->x.vals);

    phase = 331;
    pardiso(data->pt, &one, &one, &data->mtype, &phase, &data->size_matrix_mklint, A->vals, A->ptrs, A->idxs, nullptr, &one, data->iparm, &data->msglvl, data->x.vals, data->y.vals, &data->error);
    if(data->error != 0) {
        eslog::error("pardiso error %d\n", (int)data->error);
    }

    std::fill_n(data->y.vals + size_A11, size_sc, T{0});

    phase = 333;
    pardiso(data->pt, &one, &one, &data->mtype, &phase, &data->size_matrix_mklint, A->vals, A->ptrs, A->idxs, nullptr, &one, data->iparm, &data->msglvl, data->y.vals, data->x.vals, &data->error);
    if(data->error != 0) {
        eslog::error("pardiso error %d\n", (int)data->error);
    }

    std::copy_n(data->x.vals, size_A11, sol.vals);
}



template<typename T, typename I>
void schur_csx_dny_mklpardiso<T,I>::internal_solve_A11(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol)
{
    // no support for leading dimension in pardiso. matrix has to be colmajor. so I have to reallocate and convert
    Allocator_new * ator = AllocatorCPU_new::get_singleton();
    if((&rhs == &sol) && ((rhs.order != 'C') || (rhs.ld != rhs.nrows))) {
        T * tmp_mem = ator->template alloc<T>(rhs.nrows * rhs.ncols);
        MatrixDenseView_new<T> tmp;
        tmp.set_view(rhs.nrows, rhs.ncols, rhs.nrows, 'C', tmp_mem, ator);
        math::operations::convert_dnx_dny<T>::do_all(&rhs, &tmp, false);
        this->internal_solve_A11(tmp, tmp);
        math::operations::convert_dnx_dny<T>::do_all(&tmp, &sol, false);
        ator->free(tmp_mem);
        return;
    }
    if((rhs.order != 'C') || (rhs.ld != rhs.nrows)) {
        T * rhs_2_mem = ator->template alloc<T>(rhs.nrows * rhs.ncols);
        MatrixDenseView_new<T> rhs_2;
        rhs_2.set_view(rhs.nrows, rhs.ncols, rhs.nrows, 'C', rhs_2_mem, ator);
        math::operations::convert_dnx_dny<T>::do_all(&rhs, &rhs_2, false);
        this->internal_solve_A11(rhs_2, sol);
        ator->free(rhs_2_mem);
        return;
    }
    if((sol.order != 'C') || (sol.ld != sol.nrows)) {
        T * sol_2_mem = ator->template alloc<T>(sol.nrows * sol.ncols);
        MatrixDenseView_new<T> sol_2;
        sol_2.set_view(sol.nrows, sol.ncols, sol.nrows, 'C', sol_2_mem, ator);
        this->internal_solve_A11(rhs, sol_2);
        math::operations::convert_dnx_dny<T>::do_all(&sol_2, &sol, false);
        ator->free(sol_2_mem);
        return;
    }

    T * tmp_mem = ator->template alloc<T>(rhs.nrows * rhs.ncols);
    MatrixDenseView_new<T> tmp;
    tmp.set_view(rhs.nrows, rhs.ncols, rhs.nrows, 'C', tmp_mem, ator);

    MKL_INT phase;
    MKL_INT one = 1;
    MKL_INT nrhs = rhs.ncols;

    phase = 331;
    pardiso(data->pt, &one, &one, &data->mtype, &phase, &data->size_matrix_mklint, A->vals, A->ptrs, A->idxs, nullptr, &nrhs, data->iparm, &data->msglvl, rhs.vals, tmp.vals, &data->error);
    if(data->error != 0) {
        eslog::error("pardiso error %d\n", (int)data->error);
    }

    {
        MatrixDenseView_new<T> tmp_sub_bottom = tmp.get_submatrix_view(size_A11, size_matrix, 0, tmp.ncols);
        math::operations::fill_dnx<T>::do_all(&tmp_sub_bottom, T{0});
    }

    phase = 333;
    pardiso(data->pt, &one, &one, &data->mtype, &phase, &data->size_matrix_mklint, A->vals, A->ptrs, A->idxs, nullptr, &nrhs, data->iparm, &data->msglvl, tmp.vals, sol.vals, &data->error);
    if(data->error != 0) {
        eslog::error("pardiso error %d\n", (int)data->error);
    }

    ator->free(tmp_mem);
}



#define INSTANTIATE_T_I(T,I) \
template class schur_csx_dny_mklpardiso<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

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

#endif
