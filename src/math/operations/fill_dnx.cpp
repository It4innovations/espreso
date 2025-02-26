
#include "math/operations/fill_dnx.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
void fill_dnx<T>::set_matrix(MatrixDenseView_new<T> * M_)
{
    M = M_;
}



template<typename T>
void fill_dnx<T>::set_value(T val_)
{
    val = val_;
}



template<typename T>
void fill_dnx<T>::perform()
{
    if(M == nullptr) eslog::error("matrix is not set\n");
    if((M->prop.uplo == 'U' || M->prop.uplo == 'L') && M->nrows != M->ncols) eslog::error("uplo matrix must be square\n");

    size_t num_blocks = M->get_num_blocks();
    size_t block_size = M->get_block_size();
    bool move_start = ((M->uplo == 'U' && M->order == 'R') || (M->uplo == 'L' && M->order == 'C'));
    bool move_end   = ((M->uplo == 'L' && M->order == 'R') || (M->uplo == 'U' && M->order == 'C'));
    for(size_t i = 0; i < num_blocks; i++) {
        size_t start = 0;
        size_t end = block_size;
        if(move_start) start = i;
        if(move_end) end = i;
        size_t size = end - start;
        std::fill_n(M->vals + i * M->ld + start, size, val);
    }
}



template<typename T>
void fill_dnx<T>::do_all(MatrixDenseView_new<T> * M, T val)
{
    fill_dnx<T> instance;
    instance.set_matrix(M);
    instance.set_value(val);
    instance.perform();
}




}
}
}
