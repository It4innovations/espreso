
#ifndef SRC_MATH_PRIMITIVES_NEW_MATRIX_BASE_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MATRIX_BASE_NEW_H_

#include <cstddef>



namespace espreso {



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
        // TODO: symmetric, positive definite, ...
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



}



#endif /* SRC_MATH_PRIMITIVES_NEW_MATRIX_BASE_NEW_H_ */
