
#include "math/operations/fill_dnx.h"

#include "basis/utilities/stacktimer.h"



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
    stacktimer::push("fill_dnx::perform");

    if(M == nullptr) eslog::error("matrix is not set\n");
    if((M->prop.uplo == 'U' || M->prop.uplo == 'L') && M->nrows != M->ncols) eslog::error("uplo matrix must be square\n");

    size_t size_primary = M->get_size_primary();
    size_t size_secdary = M->get_size_secdary();
    bool move_start = ((M->prop.uplo == 'U' && M->order == 'R') || (M->prop.uplo == 'L' && M->order == 'C'));
    bool move_end   = ((M->prop.uplo == 'L' && M->order == 'R') || (M->prop.uplo == 'U' && M->order == 'C'));
    for(size_t i = 0; i < size_primary; i++) {
        size_t start = 0;
        size_t end = size_secdary;
        if(move_start) start = i;
        if(move_end) end = i+1;
        size_t size = end - start;
        std::fill_n(M->vals + i * M->ld + start, size, val);
    }

    stacktimer::pop();
}



template<typename T>
void fill_dnx<T>::do_all(MatrixDenseView_new<T> * M, T val)
{
    fill_dnx<T> instance;
    instance.set_matrix(M);
    instance.set_value(val);
    instance.perform();
}



#define INSTANTIATE_T(T) \
template class fill_dnx<T>;

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
