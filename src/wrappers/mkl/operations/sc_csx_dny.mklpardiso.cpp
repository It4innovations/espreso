
#ifdef HAVE_MKL

#include "wrappers/mkl/operations/sc_csx_dny.mklpardiso.h"

#include <mkl.h>

#include "math/primitives_new/allocator_new.h"
#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/vector_dense_data_new.h"
#include "math/operations/concat_csx.h"
#include "math/operations/lincomb_matrix_dnx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
struct sc_csx_dny_mklpardiso_data
{
    MatrixCsxData_new<T,I> A_whole;
    VectorDenseData_new<MKL_INT> perm;
    VectorDenseData_new<T> x;
    VectorDenseData_new<T> y;
    void * pt[64];
    MKL_INT iparm[64];
    MKL_INT mtype;
    MKL_INT msglvl;
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
sc_csx_dny_mklpardiso<T,I>::sc_csx_dny_mklpardiso() = default;



template<typename T, typename I>
sc_csx_dny_mklpardiso<T,I>::~sc_csx_dny_mklpardiso() = default;



template<typename T, typename I>
void sc_csx_dny_mklpardiso<T,I>::internal_preprocess()
{
    if(called_set_matrix == '1' && A->order != 'R') eslog::error("single input matrix needs to be rowmajor\n");

    data = std::make_unique<sc_csx_dny_mklpardiso_data<T,I>>();

    if(size_matrix > (size_t)std::numeric_limits<MKL_INT>::max()) eslog::error("size_matrix too large for MKL_INT\n");
    data->size_matrix_mklint = (MKL_INT)size_matrix;

    if(called_set_matrix == '4') {
        if(is_matrix_hermitian) {
            if(A11->prop.uplo != 'U') eslog::error("only uplo=U is supported\n");
            if(A22 != nullptr && A22->prop.uplo != 'U') eslog::error("only uplo=U is supported\n");
            if(A21 != nullptr) eslog::error("A21 should be empty\n");
            if(A11->prop.dfnt != MatrixDefinitness_new::positive_definite) eslog::error("A11 should be positive definite\n");
        }
        size_t total_nnz = 0;
        if(A11 != nullptr) total_nnz += A11->nnz;
        if(A12 != nullptr) total_nnz += A12->nnz;
        if(A21 != nullptr) total_nnz += A21->nnz;
        if(A22 != nullptr) total_nnz += A22->nnz;
        data->A_whole.set(size_matrix, size_matrix, total_nnz, 'R', AllocatorCPU_new::get_singleton());
        data->A_whole.alloc();
        if(is_matrix_hermitian) {
            data->A_whole.prop.uplo = 'U';
            data->A_whole.prop.symm = MatrixSymmetry_new::hermitian;
            data->A_whole.prop.dfnt = MatrixDefinitness_new::positive_semidefinite;
        }
        else {
            data->A_whole.prop.uplo = 'F';
            data->A_whole.prop.symm = MatrixSymmetry_new::general;
            data->A_whole.prop.dfnt = MatrixDefinitness_new::indefinite;
        }
        A = &data->A_whole;
        concat_csx<T,I>::do_all(A11, A12, A21, A22, A);
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
void sc_csx_dny_mklpardiso<T,I>::internal_perform_1()
{
    if(called_set_matrix == '4') {
        concat_csx<T,I>::do_all(A11, A12, A21, A22, A);
    }

    T * sc_tmp_vals = new T[size_sc * size_sc];
    MatrixDenseView_new<T> sc_tmp;
    sc_tmp.set_view(size_sc, size_sc, size_sc, sc->order, sc_tmp_vals, AllocatorGPU_new::get_singleton());

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
void sc_csx_dny_mklpardiso<T,I>::internal_perform_2()
{
}



template<typename T, typename I>
void sc_csx_dny_mklpardiso<T,I>::internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
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



#define INSTANTIATE_T_I(T,I) \
template class sc_csx_dny_mklpardiso<T,I>;

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
