
#ifndef SRC_MATH_PRIMITIVES_NEW_MATRIX_BASE_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MATRIX_BASE_NEW_H_



template<typename T>
struct MatrixBase_new
{
    size_t nrows = 0;
    size_t ncols = 0;
    char diag = '_'; // Unit, Nonunit
    char uplo = '_'; // Upper, Lower

    MatrixBase_new() {}
    virtual ~MatrixBase_new() {}

    MatrixBase_new(const MatrixBase_new &) = default;
    MatrixBase_new(MatrixBase_new &&) = default;
    MatrixBase_new & operator=(const MatrixBase_new &) = default;
    MatrixBase_new & operator=(MatrixBase_new &&) = default;
};



char change_uplo(char uplo)
{
    switch(uplo)
    {
        case 'U': return 'L';
        case 'L': return 'U';
        case 'F': return 'F';
        default: return '_';
    }
}

char change_order(char order)
{
    switch(order)
    {
        case 'R': return 'C';
        case 'C': return 'R';
        default: return '_';
    }
}



#endif /* SRC_MATH_PRIMITIVES_NEW_MATRIX_BASE_NEW_H_ */
