
#include "math/operations/concat_csx.h"

#include "math/primitives_new/allocator_new.h"
#include "math/primitives_new/matrix_csx_data_new.h"
#include "math/operations/convert_csx_csy.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void concat_csx<T,I>::set_matrices_src(MatrixCsxView_new<T,I> * A11_, MatrixCsxView_new<T,I> * A12_, MatrixCsxView_new<T,I> * A21_, MatrixCsxView_new<T,I> * A22_)
{
    if(A11 != nullptr || A12 != nullptr || A21 != nullptr || A22 != nullptr) eslog::error("source matrices are already set\n");

    A11 = A11_;
    A12 = A12_;
    A21 = A21_;
    A22 = A22_;
}



template<typename T, typename I>
void concat_csx<T,I>::set_matrix_dst(MatrixCsxView_new<T,I> * A_)
{
    if(A != nullptr) eslog::error("matrix A is already set\n");

    A = A_;
}



template<typename T, typename I>
void concat_csx<T,I>::perform()
{
    stacktimer::push("concat_csx::perform");

    if(A == nullptr) eslog::error("destination matrix is not set\n");
    if(A11 != nullptr && !A11->ator->is_data_accessible_cpu()) eslog::error("matrix A11 must be cpu-accessible\n");
    if(A12 != nullptr && !A12->ator->is_data_accessible_cpu()) eslog::error("matrix A12 must be cpu-accessible\n");
    if(A21 != nullptr && !A21->ator->is_data_accessible_cpu()) eslog::error("matrix A21 must be cpu-accessible\n");
    if(A22 != nullptr && !A22->ator->is_data_accessible_cpu()) eslog::error("matrix A22 must be cpu-accessible\n");
    if(!A->ator->is_data_accessible_cpu()) eslog::error("matrix A must be cpu-accessible\n");
    size_t total_nnz = 0;
    if(A11 != nullptr) total_nnz += A11->nnz;
    if(A12 != nullptr) total_nnz += A12->nnz;
    if(A21 != nullptr) total_nnz += A21->nnz;
    if(A22 != nullptr) total_nnz += A22->nnz;
    if(total_nnz != A->nnz) eslog::error("wrong nnz of output matrix\n");
    size_t nrows_top = 0;
    size_t nrows_bot = 0;
    size_t ncols_left = 0;
    size_t ncols_right = 0;
    if(A11 != nullptr) {
        nrows_top = A11->nrows;
        nrows_bot = A->nrows - A11->nrows;
        ncols_left = A11->ncols;
        ncols_right = A->ncols - A11->ncols;
    }
    else if(A12 != nullptr) {
        nrows_top = A12->nrows;
        nrows_bot = A->nrows - A12->nrows;
        ncols_left = A->ncols - A12->ncols;
        ncols_right = A12->ncols;
    }
    else if(A21 != nullptr) {
        nrows_top = A->nrows - A21->nrows;
        nrows_bot = A21->nrows;
        ncols_left = A21->ncols;
        ncols_right = A->ncols - A21->ncols;
    }
    else if(A22 != nullptr) {
        nrows_top = A->nrows - A22->nrows;
        nrows_bot = A22->nrows;
        ncols_left = A->ncols - A22->ncols;
        ncols_right = A22->ncols;
    }
    else {
        eslog::error("at least one input matrix must be set\n");
    }
    if(A11 != nullptr && (A11->nrows != nrows_top || A11->ncols != ncols_left)) eslog::error("wrong size of A11\n");
    if(A12 != nullptr && (A12->nrows != nrows_top || A12->ncols != ncols_right)) eslog::error("wrong size of A12\n");
    if(A21 != nullptr && (A21->nrows != nrows_bot || A21->ncols != ncols_left)) eslog::error("wrong size of A21\n");
    if(A22 != nullptr && (A22->nrows != nrows_bot || A22->ncols != ncols_right)) eslog::error("wrong size of A22\n");

    MatrixCsxData_new<T,I> A11_tmp;
    MatrixCsxData_new<T,I> A12_tmp;
    MatrixCsxData_new<T,I> A21_tmp;
    MatrixCsxData_new<T,I> A22_tmp;

    MatrixCsxView_new<T,I> * A11_to_use = A11;
    MatrixCsxView_new<T,I> * A12_to_use = A12;
    MatrixCsxView_new<T,I> * A21_to_use = A21;
    MatrixCsxView_new<T,I> * A22_to_use = A22;

    if(A11 != nullptr && A11->order != A->order) {
        A11_tmp.set(A11->nrows, A11->ncols, A11->nnz, A->order, AllocatorCPU_new::get_singleton());
        A11_tmp.alloc();
        convert_csx_csy<T,I>::do_all(A11, &A11_tmp);
        A11_to_use = &A11_tmp;
    }
    if(A12 != nullptr && A12->order != A->order) {
        A12_tmp.set(A12->nrows, A12->ncols, A12->nnz, A->order, AllocatorCPU_new::get_singleton());
        A12_tmp.alloc();
        convert_csx_csy<T,I>::do_all(A12, &A12_tmp);
        A12_to_use = &A12_tmp;
    }
    if(A21 != nullptr && A21->order != A->order) {
        A21_tmp.set(A21->nrows, A21->ncols, A21->nnz, A->order, AllocatorCPU_new::get_singleton());
        A21_tmp.alloc();
        convert_csx_csy<T,I>::do_all(A21, &A21_tmp);
        A21_to_use = &A21_tmp;
    }
    if(A22 != nullptr && A22->order != A->order) {
        A22_tmp.set(A22->nrows, A22->ncols, A22->nnz, A->order, AllocatorCPU_new::get_singleton());
        A22_tmp.alloc();
        convert_csx_csy<T,I>::do_all(A22, &A22_tmp);
        A22_to_use = &A22_tmp;
    }

    I * const A_ptrs = A->ptrs;
    I * const A_idxs = A->idxs;
    T * const A_vals = A->vals;

    size_t midpoint_primary = ((A->order == 'R') ? nrows_top : ncols_left);
    size_t midpoint_secdary = ((A->order == 'R') ? ncols_left : nrows_top);

    size_t curr_idx = 0;

    auto incorporate_matrix_part = [&](MatrixCsxView_new<T,I> const * const M, size_t ipd, size_t offset_primary, size_t offset_secdary){
        if(M != nullptr) {
            size_t ips = ipd - offset_primary;
            I start = M->ptrs[ips];
            I end = M->ptrs[ips+1];
            for(I i = start; i < end; i++) {
                I iss = M->idxs[i];
                I isd = iss + offset_secdary;
                T val = M->vals[i];
                A_idxs[curr_idx] = isd;
                A_vals[curr_idx] = val;
                curr_idx++;
            }
        }
    };

    for(size_t ipd = 0; ipd < A->get_size_primary(); ipd++) {
        A_ptrs[ipd] = curr_idx;
        if(ipd < midpoint_primary) {
            incorporate_matrix_part(A11_to_use, ipd, 0, 0);
            if(A->order == 'R') incorporate_matrix_part(A12_to_use, ipd, 0, midpoint_secdary);
            if(A->order == 'C') incorporate_matrix_part(A21_to_use, ipd, 0, midpoint_secdary);
        }
        else {
            if(A->order == 'R') incorporate_matrix_part(A21_to_use, ipd, midpoint_primary, 0);
            if(A->order == 'C') incorporate_matrix_part(A12_to_use, ipd, midpoint_primary, 0);
            incorporate_matrix_part(A22_to_use, ipd, midpoint_primary, midpoint_secdary);
        }
    }
    A_ptrs[A->get_size_primary()] = curr_idx;

    stacktimer::pop();
}



template<typename T, typename I>
void concat_csx<T,I>::do_all(MatrixCsxView_new<T,I> * A11, MatrixCsxView_new<T,I> * A12, MatrixCsxView_new<T,I> * A21, MatrixCsxView_new<T,I> * A22, MatrixCsxView_new<T,I> * A)
{
    concat_csx<T,I> instance;
    instance.set_matrices_src(A11, A12, A21, A22);
    instance.set_matrix_dst(A);
    instance.perform();
}



#define INSTANTIATE_T_I(T,I) \
template class concat_csx<T,I>;

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
