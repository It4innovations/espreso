
#ifndef SRC_WRAPPERS_MKL_PARDISO_COMMON_H
#define SRC_WRAPPERS_MKL_PARDISO_COMMON_H

#include <mkl.h>

#include "math/primitives_new/matrix_base_new.h"



namespace espreso {



template<typename T>
inline int get_pardiso_matrix_type(const MatrixBase_new::matrix_properties & prop)
{
    if constexpr(utils::is_complex<T>()) {
        if(prop.symm == MatrixSymmetry_new::hermitian) {
            if(prop.dfnt == MatrixDefinitness_new::positive_definite) {
                return 4; // complex and Hermitian positive definite
            }
            else {
                return -4; // complex and Hermitian indefinite
                // we have other MatrixDefinitness_new values, but pardiso does not have them, so we will treat it as indefinite
            }
        }
        if(prop.symm == MatrixSymmetry_new::symmetric) {
            return 6; // complex and symmetric
        }
        if(prop.symm == MatrixSymmetry_new::structurally_symmetric) {
            return 3; // complex and structurally symmetric
        }
        if(prop.symm == MatrixSymmetry_new::general) {
            return 13; // complex and nonsymmetric
        }
    }
    if constexpr(utils::is_real<T>()) {
        if(is_symmetric<T>(prop.symm)) {
            if(prop.dfnt == MatrixDefinitness_new::positive_definite) {
                return 2; // real and symmetric positive definite
            }
            else {
                return -2; // real and symmetric indefinite
                // we have other MatrixDefinitness_new values, but pardiso does not have them, so we will treat it as indefinite
            }
        }
        if(prop.symm == MatrixSymmetry_new::structurally_symmetric) {
            return 1; // real and structurally symmetric
        }
        if(prop.symm == MatrixSymmetry_new::general) {
            return 11; // real and nonsymmetric
        }
    }

    eslog::error("weird matrix\n");
}



inline const char * get_pardiso_error_string(MKL_INT error)
{
    switch (error) {
    case   0: return "success";
    case  -1: return "input inconsistent";
    case  -2: return "not enough memory";
    case  -3: return "reordering problem";
    case  -4: return "zero pivot, numerical factorization or iterative refinement problem";
    case  -5: return "unclassified (internal) error";
    case  -6: return "reordering failed";
    case  -7: return "diagonal matrix is singular";
    case  -8: return "32-bit integer overflow problem";
    case  -9: return "not enough memory for OOC";
    case -10: return "error opening OOC files";
    case -11: return "read/write error with OOC files";
    case -12: return "(pardiso_64 only) pardiso_64 called from 32-bit library";
    case -13: return "interrupted by the (user-defined) mkl_progress function";
    case -15: return "internal error which can appear for iparm[23]=10 and iparm[12]=1. Try switch matching off (set iparm[12]=0 and rerun";
    }
    return "unspecified pardiso error";
}



}



#endif /* SRC_WRAPPERS_MKL_PARDISO_COMMON_H */
