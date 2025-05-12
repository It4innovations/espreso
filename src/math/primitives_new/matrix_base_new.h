
#ifndef SRC_MATH_PRIMITIVES_NEW_MATRIX_BASE_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MATRIX_BASE_NEW_H_

#include <cstddef>
#include <cstdint>

#include "basis/utilities/utils.h"
#include "math/primitives/matrix_info.h"



namespace espreso {



enum struct MatrixUplo_new : uint8_t
{
    unset = 0,
    full,
    upper,
    lower
};

enum struct MatrixDiag_new : uint8_t
{
    unset = 0,
    nonunit,
    unit
};

enum struct MatrixSymmetry_new : uint8_t
{
    unset = 0,
    general,
    structurally_symmetric,
    symmetric,
    hermitian
};

enum struct MatrixDefinitness_new : uint8_t
{
    unset = 0,
    indefinite,
    positive_definite,
    positive_semidefinite,
    negative_definite,
    negative_semidefinite
};



class MatrixBase_new
{
public: // the user promises not to modify these values (I don't want to implement getters everywhere)
    size_t nrows = 0;
    size_t ncols = 0;
public: // ok to modify
    struct matrix_properties
    {
        char uplo = '_'; // Upper, Lower
        char diag = '_'; // Unit, Nonunit
        MatrixSymmetry_new symm = MatrixSymmetry_new::unset;
        MatrixDefinitness_new dfnt = MatrixDefinitness_new::unset;
    } prop;
public:
    MatrixBase_new() = default;
    MatrixBase_new(const MatrixBase_new &) = default;
    MatrixBase_new(MatrixBase_new &&) = default;
    MatrixBase_new & operator=(const MatrixBase_new &) = default;
    MatrixBase_new & operator=(MatrixBase_new &&) = default;
    virtual ~MatrixBase_new() = default;
};



inline char change_uplo(char uplo)
{
    switch(uplo)
    {
        case 'U': return 'L';
        case 'L': return 'U';
        case 'F': return 'F';
        default: return '_';
    }
}

inline char change_order(char order)
{
    switch(order)
    {
        case 'R': return 'C';
        case 'C': return 'R';
        default: return '_';
    }
}

template<typename T>
inline bool is_symmetric(MatrixSymmetry_new symm)
{
    // is it true, that A == At ?
    if(symm == MatrixSymmetry_new::symmetric) return true;
    if constexpr(utils::is_real<T>()) if(symm == MatrixSymmetry_new::hermitian) return true;
    return false;
}

template<typename T>
inline bool is_hermitian(MatrixSymmetry_new symm)
{
    // is it true, that A == A* ?
    if(symm == MatrixSymmetry_new::hermitian) return true;
    if constexpr(utils::is_real<T>()) if(symm == MatrixSymmetry_new::symmetric) return true;
    return false;
}

template<typename T>
inline bool is_symmetry_equal(MatrixSymmetry_new s1, MatrixSymmetry_new s2)
{
    if constexpr(utils::is_real<T>()) return is_symmetric<T>(s1) == is_symmetric<T>(s2);
    return s1 == s2;
}

inline Matrix_Shape get_old_matrix_shape(const MatrixBase_new & A)
{
    if(A.prop.uplo == 'U') return Matrix_Shape::UPPER;
    else if(A.prop.uplo == 'L') return Matrix_Shape::LOWER;
    else return Matrix_Shape::FULL;
}

template<typename T>
inline Matrix_Type get_old_matrix_type(const MatrixBase_new & A)
{
    if constexpr(utils::is_complex<T>()) {
        if(A.prop.symm == MatrixSymmetry_new::hermitian && A.prop.dfnt == MatrixDefinitness_new::positive_definite) return Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE;
        if(A.prop.symm == MatrixSymmetry_new::hermitian && A.prop.dfnt == MatrixDefinitness_new::indefinite) return Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE;
        if(A.prop.symm == MatrixSymmetry_new::hermitian) return Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE;  // not really true, but no other matrix_type options
        if(A.prop.symm == MatrixSymmetry_new::symmetric) return Matrix_Type::COMPLEX_SYMMETRIC;
        if(A.prop.symm == MatrixSymmetry_new::structurally_symmetric) return Matrix_Type::COMPLEX_STRUCTURALLY_SYMMETRIC;
        if(A.prop.symm == MatrixSymmetry_new::general) return Matrix_Type::COMPLEX_NONSYMMETRIC;
    }
    if constexpr(utils::is_real<T>()) {
        if(is_symmetric<T>(A.prop.symm) && A.prop.dfnt == MatrixDefinitness_new::positive_definite) return Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
        if(is_symmetric<T>(A.prop.symm) && A.prop.dfnt == MatrixDefinitness_new::indefinite) return Matrix_Type::REAL_SYMMETRIC_INDEFINITE;
        if(is_symmetric<T>(A.prop.symm)) return Matrix_Type::REAL_SYMMETRIC_INDEFINITE; // not really true, but no other matrix_type options
        if(A.prop.symm == MatrixSymmetry_new::structurally_symmetric) return Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
        if(A.prop.symm == MatrixSymmetry_new::general) return Matrix_Type::REAL_NONSYMMETRIC;
    }

    return Matrix_Type::UNSET_INVALID_NONE;
}

inline char get_new_matrix_uplo(Matrix_Shape ms)
{
    if(ms == Matrix_Shape::UPPER) return 'U';
    if(ms == Matrix_Shape::LOWER) return 'L';
    if(ms == Matrix_Shape::FULL) return 'F';
    return '_';
}

inline MatrixSymmetry_new get_new_matrix_symmetry(Matrix_Type mt)
{
    switch(mt) {
        case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
        case Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE:
            return MatrixSymmetry_new::hermitian;
        case Matrix_Type::COMPLEX_SYMMETRIC:
        case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
        case Matrix_Type::REAL_SYMMETRIC_INDEFINITE:
            return MatrixSymmetry_new::symmetric;
        case Matrix_Type::COMPLEX_STRUCTURALLY_SYMMETRIC:
        case Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC:
            return MatrixSymmetry_new::structurally_symmetric;
        case Matrix_Type::COMPLEX_NONSYMMETRIC:
        case Matrix_Type::REAL_NONSYMMETRIC:
            return MatrixSymmetry_new::general;
        default:
            return MatrixSymmetry_new::unset;
    }
}

inline MatrixDefinitness_new get_new_matrix_definitness(Matrix_Type mt)
{
    switch(mt) {
        case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
        case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
            return MatrixDefinitness_new::positive_definite;
        case Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE:
        case Matrix_Type::REAL_SYMMETRIC_INDEFINITE:
            return MatrixDefinitness_new::indefinite;
        case Matrix_Type::COMPLEX_SYMMETRIC:
        case Matrix_Type::COMPLEX_STRUCTURALLY_SYMMETRIC:
        case Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC:
        case Matrix_Type::COMPLEX_NONSYMMETRIC:
        case Matrix_Type::REAL_NONSYMMETRIC:
        default:
            return MatrixDefinitness_new::unset;
    }
}



}



#endif /* SRC_MATH_PRIMITIVES_NEW_MATRIX_BASE_NEW_H_ */
