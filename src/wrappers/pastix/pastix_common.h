
#ifndef SRC_WRAPPERS_PASTIX_PASTIX_COMMON_H
#define SRC_WRAPPERS_PASTIX_PASTIX_COMMON_H



#include <pastix.h>

#include "math/primitives_new/matrix_base_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T>
spm_mtxtype_t symm_to_pastix(MatrixSymmetry_new symm)
{
    // I have to set hermitian only if I use complex numbers, for real I have to use symmetric
    if(utils::is_complex<T>() && is_hermitian<T>(symm)) {
        return SpmHermitian;
    }
    else if(is_symmetric<T>(symm)) {
        return SpmSymmetric;
    }
    else {
        return SpmGeneral;
    }
}



template<typename T>
spm_coeftype_t type_to_pastix()
{
    if constexpr(std::is_same_v<T,float>) return SpmFloat;
    if constexpr(std::is_same_v<T,double>) return SpmDouble;
    if constexpr(std::is_same_v<T,std::complex<float>>) return SpmComplex32;
    if constexpr(std::is_same_v<T,std::complex<double>>) return SpmComplex64;
}



spm_fmttype_e order_to_pastix(char order)
{
    if(order == 'R') return SpmCSR;
    if(order == 'C') return SpmCSC;
    eslog::error("wrong order\n");
}



}
}
}



#endif /* SRC_WRAPPERS_PASTIX_PASTIX_COMMON_H */
