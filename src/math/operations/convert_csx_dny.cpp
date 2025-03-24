
#include "math/operations/convert_csx_dny.h"

#include "math/operations/fill_dnx.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void convert_csx_dny<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    M_src = M_src_;
}



template<typename T, typename I>
void convert_csx_dny<T,I>::set_matrix_dst(MatrixDenseView_new<T> * M_dst_)
{
    M_dst = M_dst_;
}



template<typename T, typename I>
void convert_csx_dny<T,I>::perform_zerofill()
{
    stacktimer::push("convert_csx_dny::perform_zerofill");

    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix dimensions don't match\n");
    if(M_src->prop.uplo != M_dst->prop.uplo) eslog::error("uplo of matrices does not match\n");
    if((M_dst->prop.uplo == 'L' || M_dst->prop.uplo == 'U') && M_dst->nrows != M_dst->ncols) eslog::error("upper of lower matrix must be square\n");

    fill_dnx<T>::do_all(M_dst, T{0});

    stacktimer::pop();
}



template<typename T, typename I>
void convert_csx_dny<T,I>::perform_copyvals()
{
    stacktimer::push("convert_csx_dny::perform_copyvals");

    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix dimensions don't match\n");
    if(M_src->prop.uplo != M_dst->prop.uplo) eslog::error("uplo of matrices does not match\n");
    if((M_dst->prop.uplo == 'L' || M_dst->prop.uplo == 'U') && M_dst->nrows != M_dst->ncols) eslog::error("upper of lower matrix must be square\n");

    size_t primary_size = M_src->get_size_primary();
    size_t dstld = M_dst->ld;
    I * srcptrs = M_src->ptrs;
    I * srcidxs = M_src->idxs;
    T * srcvals = M_src->vals;
    T * dstvals = M_dst->vals;

    if(M_src->order == M_dst->order)
    {
        for(size_t r = 0; r < primary_size; r++)
        {
            I start = srcptrs[r];
            I end = srcptrs[r+1];
            T * row = dstvals + r * dstld;
            for(I i = start; i < end; i++)
            {
                I c = srcidxs[i];
                T v = srcvals[i];
                row[c] = v;
            }
        }
    }
    else
    {
        for(size_t r = 0; r < primary_size; r++)
        {
            I start = srcptrs[r];
            I end = srcptrs[r+1];
            for(I i = start; i < end; i++)
            {
                I c = srcidxs[i];
                T v = srcvals[i];
                dstvals[r + dstld * c] = v;
            }
        }
    }

    stacktimer::pop();
}



template<typename T, typename I>
void convert_csx_dny<T,I>::perform_all()
{
    perform_zerofill();
    perform_copyvals();
}



template<typename T, typename I>
void convert_csx_dny<T,I>::do_all(MatrixCsxView_new<T,I> * M_src, MatrixDenseView_new<T> * M_dst)
{
    convert_csx_dny<T,I> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.perform_zerofill();
    instance.perform_copyvals();
}



#define INSTANTIATE_T_I(T,I) \
template class convert_csx_dny<T,I>;

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
