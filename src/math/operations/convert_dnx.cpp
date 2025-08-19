
#include "math/operations/convert_dnx.h"

#include "math/primitives_new.h"
#include "math/primitives_new/allocator_new.h"
#include "math/operations/copy_dnx.h"
#include "math/operations/convert_dnx_dny.h"
#include "math/operations/complete_dnx_dnx.h"
#include "basis/utilities/stacktimer.h"


namespace espreso {
namespace math {
namespace operations {



template<typename T>
void convert_dnx<T>::set_matrix_src(MatrixDenseView_new<T> * M_src_)
{
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");

    M_src = M_src_;
}



template<typename T>
void convert_dnx<T>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("matrix M_dst is already set\n");

    M_dst = M_dst_;
}



template<typename T>
void convert_dnx<T>::perform()
{
    stacktimer::push("convert_dnx_dny::perform");

    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_src == M_dst) eslog::error("in-place makes no sense here\n");
    if(!M_src->ator->is_data_accessible_cpu()) eslog::error("source matrix must be cpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_cpu()) eslog::error("destination matrix must be cpu-accessible\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");
    if(M_src->prop.symm != M_dst->prop.symm) eslog::error("matrix symemtry does not match\n");
    MatrixSymmetry_new symm = M_src->prop.symm;
    if(is_structurally_symmetric(symm) && M_src->nrows != M_src->ncols) eslog::error("symmetric matrix must be square\n");

    const char src_order = M_src->order;
    const char dst_order = M_dst->order;
    const char src_uplo = M_src->prop.uplo;
    const char dst_uplo = M_dst->prop.uplo;

    if(is_structurally_symmetric(symm)) {
        bool herm = (utils::is_complex<T>() && symm == MatrixSymmetry_new::hermitian);
        if(src_order == dst_order) {
            // same order
            if(is_uplo_equal(src_uplo, dst_uplo)) {
                // same order and same uplo, just copy
                // L->L, U->U, F->F
                copy_dnx<T>::do_all(M_src, M_dst, false);
            }
            else if(is_uplo(src_uplo) && !is_uplo(dst_uplo)) {
                // L->F, U->F
                // we need to complete the matrix
                complete_dnx_dnx<T>::do_all(M_src, M_dst, herm);
            }
            else if(!is_uplo(src_uplo) && is_uplo(dst_uplo)) {
                // F->L, F->U
                // extract just a single triangle
                MatrixDenseView_new<T> M_tmp = *M_src;
                M_tmp.prop.uplo = dst_uplo;
                copy_dnx<T>::do_all(&M_tmp, M_dst);
            }
            else {
                // L->U, U->L
                // M_tmp - uplo same as src, order different than src=dst
                // M_tmp_rt - uplo same as dst, order same as dst
                MatrixDenseData_new<T> M_tmp;
                M_tmp.set(M_src->nrows, M_src->ncols, change_order(src_order), AllocatorCPU_new::get_singleton());
                M_tmp.prop.uplo = src_uplo;
                M_tmp.alloc();
                convert_dnx_dny<T>::do_all(M_src, &M_tmp, false);
                MatrixDenseView_new<T> M_tmp_rt = M_tmp.get_transposed_reordered_view();
                copy_dnx<T>::do_all(&M_tmp_rt, M_dst, herm);
            }
        }
        else {
            // different order
            if(is_uplo_equal(src_uplo, dst_uplo)) {
                // L->L, U->U, F->F
                // need to only reorder, but convert_dnx_dny ignores uplo, so I might neet to go through a tmp matrix
                if(is_uplo(src_uplo)) {
                    // L->L, U->U
                    // M_tmp - uplo from src=dst, order from dst
                    MatrixDenseData_new<T> M_tmp;
                    M_tmp.set(M_src->nrows, M_src->ncols, dst_order, AllocatorCPU_new::get_singleton());
                    M_tmp.prop.uplo = src_uplo;
                    M_tmp.alloc();
                    convert_dnx_dny<T>::do_all(M_src, &M_tmp, false);
                    copy_dnx<T>::do_all(&M_tmp, M_dst, false);
                }
                else {
                    // F->F
                    convert_dnx_dny<T>::do_all(M_src, M_dst, false);
                }
            }
            else if(is_uplo(src_uplo) && !is_uplo(dst_uplo)) {
                // L->F, U->F
                // first make a same-order view of src, then just complete
                MatrixDenseView_new<T> M_src_rt = M_src->get_transposed_reordered_view();
                complete_dnx_dnx<T>::do_all(&M_src_rt, M_dst, herm);
                if constexpr(utils::is_complex<T>()) if(herm) {
                    size_t size_dst_primary = M_dst->get_size_primary();
                    size_t size_dst_secdary = M_dst->get_size_secdary();
                    for(size_t ip = 0; ip < size_dst_primary; ip++) {
                        T * prim = M_dst->vals + ip * M_dst->ld;
                        for(size_t is = 0; is < size_dst_secdary; is++) {
                            prim[is] = std::conj(is);
                        }
                    }
                }
            }
            else if(!is_uplo(src_uplo) && is_uplo(dst_uplo)) {
                // F->L, F->U
                // make a same-order view of src, then just copy the triangle
                MatrixDenseView_new<T> M_src_rt = M_src->get_transposed_reordered_view();
                M_src_rt.prop.uplo = M_dst->prop.uplo;
                copy_dnx<T>::do_all(&M_src_rt, M_dst, herm);
            }
            else {
                // L->U, U->L
                // just a different view of the same matrix, and just copy
                MatrixDenseView_new<T> M_src_rt = M_src->get_transposed_reordered_view();
                copy_dnx<T>::do_all(&M_src_rt, M_dst, herm);
            }
        }
    }
    else {
        // unsymmetric, full matrix is present in both src and dst matrices
        // I only care about order, uplo does not play a role here
        // so re-order is enough, convert_dnx_dny can do this
        convert_dnx_dny<T>::do_all(M_src, M_dst, false);
    }

    stacktimer::pop();
}



template<typename T>
void convert_dnx<T>::do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst)
{
    convert_dnx<T> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.perform();
}



#define INSTANTIATE_T(T) \
template class convert_dnx<T>;

    #define INSTANTIATE \
    /* INSTANTIATE_T(float) */ \
    INSTANTIATE_T(double) \
    /* INSTANTIATE_T(std::complex<float>) */ \
    INSTANTIATE_T(std::complex<double>)

        INSTANTIATE

    #undef INSTANTIATE
#undef INSTANTIATE_T



}
}
}
