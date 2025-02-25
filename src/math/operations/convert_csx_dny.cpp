
#include "math/operations/convert_csx_dny.h"



template<typename T, typename I>
void convert_csx_dny<T,I>::set_matrix_src(MatrixCsxView_new<T,I> & M_src_)
{
    M_src = M_src_;
}



template<typename T, typename I>
void convert_csx_dny<T,I>::set_matrix_dst(const MatrixDenseView_new<T> * M_dst_)
{
    M_dst = M_dst_;
}



template<typename T, typename I>
void convert_csx_dny<T,I>::perform_zerofill()
{
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix dimensions don't match\n");
    if(M_src->prop.uplo != M_dst->prop.uplo) eslog::error("uplo of matrices does not match\n");
    if((M_dst->uplo == 'L' || M_dst->uplo == 'U') && M_dst->nrows != M_dst->ncols) eslog::error("upper of lower matrix must be square\n");

    size_t num_blocks = M_dst->get_num_blocks();
    size_t block_size = M_dst->get_block_size();
    bool move_start = ((M_dst->uplo == 'U' && M_dst->order == 'R') || (M_dst->uplo == 'L' && M_dst->order == 'C'));
    bool move_end   = ((M_dst->uplo == 'L' && M_dst->order == 'R') || (M_dst->uplo == 'U' && M_dst->order == 'C'));
    for(size_t i = 0; i < num_blocks; i++) {
        size_t start = 0;
        size_t end = block_size;
        if(move_start) start = i;
        if(move_end) end = i;
        size_t size = end - start;
        std::fill_n(M_dst->vals + i * M_dst->ld + start, size, T{0});
    }
}



template<typename T, typename I>
void convert_csx_dny<T,I>::perform_copyvals()
{
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix dimensions don't match\n");
    if(M_src->prop.uplo != M_dst->prop.uplo) eslog::error("uplo of matrices does not match\n");
    if((M_dst->uplo == 'L' || M_dst->uplo == 'U') && M_dst->nrows != M_dst->ncols) eslog::error("upper of lower matrix must be square\n");

    size_t primary_size = M_src->get_primary_size();
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
