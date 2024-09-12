
#ifndef SRC_MATH2_PRIMITIVES_MATRIX_INFO_H_
#define SRC_MATH2_PRIMITIVES_MATRIX_INFO_H_

#include "esinfo/eslog.h"

namespace espreso {

struct Indexing {
    const static int CSR;
    const static int IJV;
};

enum struct Matrix_Type: int {
    UNSET_INVALID_NONE,

    REAL_SYMMETRIC_POSITIVE_DEFINITE,
    REAL_SYMMETRIC_INDEFINITE,
    REAL_STRUCTURALLY_SYMMETRIC,
    REAL_NONSYMMETRIC,

    COMPLEX_HERMITIAN_POSITIVE_DEFINITE,
    COMPLEX_HERMITIAN_INDEFINITE,
    COMPLEX_SYMMETRIC,
    COMPLEX_STRUCTURALLY_SYMMETRIC,
    COMPLEX_NONSYMMETRIC,
};

enum struct Matrix_Shape: int {
    LOWER,
    UPPER,
    FULL
};

enum struct Matrix_Symmetry: int {
    NONE,                    // no symmetry
    STRUCTURALLY_SYMMETRIC,    // symmetric distribution of nonzeros, but unsymmetric values
    SYMMETRIC,                // A=At
    HERMITIAN                // A=A*
};

static inline Matrix_Symmetry getSymmetry(Matrix_Type mt) {
    switch(mt) {
    case Matrix_Type::UNSET_INVALID_NONE:
        eslog::error("Invalid/unset matrix type\n");
    case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
    case Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE:
    case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:   // real symmetric is actually also hermitian
    case Matrix_Type::REAL_SYMMETRIC_INDEFINITE:
        return Matrix_Symmetry::HERMITIAN;
    case Matrix_Type::COMPLEX_SYMMETRIC:
        return Matrix_Symmetry::SYMMETRIC;
    case Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC:
    case Matrix_Type::COMPLEX_STRUCTURALLY_SYMMETRIC:
        return Matrix_Symmetry::STRUCTURALLY_SYMMETRIC;
    case Matrix_Type::REAL_NONSYMMETRIC:
    case Matrix_Type::COMPLEX_NONSYMMETRIC:
        return Matrix_Symmetry::NONE;
    default:
        return Matrix_Symmetry::NONE;
    }
}

}



#endif /* SRC_MATH2_PRIMITIVES_MATRIX_INFO_H_ */
