
#include "math/operations/copy_dnx.h"

#include "basis/utilities/utils.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void copy_dnx<T>::set_matrix_src(MatrixDenseView_new<T> * M_src_)
{
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");

    M_src = M_src_;
}



template<typename T>
void copy_dnx<T>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("matrix M_dst is already set\n");

    M_dst = M_dst_;
}



template<typename T>
void copy_dnx<T>::set_conj(bool do_conj_)
{
    do_conj = do_conj_;
}



template<typename T>
void copy_dnx<T>::perform()
{
    stacktimer::push("copy_dnx::perform");

    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!M_src->ator->is_data_accessible_cpu()) eslog::error("source matrix must be cpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_cpu()) eslog::error("destination matrix must be cpu-accessible\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");
    if(M_src->order != M_dst->order) eslog::error("matrix order does not match\n");
    if(!is_uplo_equal(M_src->prop.uplo, M_dst->prop.uplo)) eslog::error("matrix uplo does not match\n");

    size_t size_primary = M_dst->get_size_primary();
    size_t size_secdary = M_dst->get_size_secdary();
    bool move_start = ((M_dst->prop.uplo == 'U' && M_dst->order == 'R') || (M_dst->prop.uplo == 'L' && M_dst->order == 'C'));
    bool move_end   = ((M_dst->prop.uplo == 'L' && M_dst->order == 'R') || (M_dst->prop.uplo == 'U' && M_dst->order == 'C'));
    for(size_t i = 0; i < size_primary; i++) {
        size_t start = 0;
        size_t end = size_secdary;
        if(move_start) start = i;
        if(move_end) end = i+1;
        size_t size = end - start;
        if constexpr (utils::is_complex<T>()) {
            if(do_conj) {
                T * sub_src = M_src->vals + i * M_src->ld + start;
                T * sub_dst = M_dst->vals + i * M_dst->ld + start;
                for(size_t j = 0; j < size; j++) {
                    sub_dst[j] = std::conj(sub_src[j]);
                }
            }
            else {
                std::copy_n(M_src->vals + i * M_src->ld + start, size, M_dst->vals + i * M_dst->ld + start);
            }
        }
        else {
            std::copy_n(M_src->vals + i * M_src->ld + start, size, M_dst->vals + i * M_dst->ld + start);
        }
    }

    stacktimer::pop();
}



template<typename T>
void copy_dnx<T>::do_all(MatrixDenseView_new<T> * M_src, MatrixDenseView_new<T> * M_dst, bool do_conj)
{
    copy_dnx<T> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.set_conj(do_conj);
    instance.perform();
}



#define INSTANTIATE_T(T) \
template class copy_dnx<T>;

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
