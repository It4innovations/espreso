
#include "esinfo/eslog.h"

#include <algorithm>
#include <complex>
#include <type_traits>
#include <vector>

namespace espreso {

template<Matrix_Shape Shape, typename I>
static I denseMatCalcIndexRowMajor(I r, I c, I ncols)
{
    if constexpr(Shape == Matrix_Shape::LOWER) { return r * ncols + (r * (r - 1) / 2)     + c; }
    if constexpr(Shape == Matrix_Shape::FULL)  { return r * ncols                         + c; }
    if constexpr(Shape == Matrix_Shape::UPPER) { return r * ncols - (r * (r - 1) / 2) - r + c; }
    return 0;
}

template <typename T, typename I, bool DoTrans, bool DoConj, Matrix_Symmetry OutSymmetry, Matrix_Shape InShape, Matrix_Shape OutShape>
static void _submatrix(const Matrix_CSR<T, I> &input, Matrix_Dense<T, I> &output, I start_row, I end_row, I start_col, I end_col)
{
    // for upper/lower to full matrices, copies only the part that is present in the input matrix

    if constexpr(InShape == Matrix_Shape::FULL && OutShape != Matrix_Shape::FULL) {
        eslog::error("Submatrix: invalid combination of parameters - if input matrix is full, output matrix also has to be full\n");
        return;
    }
    if constexpr((OutSymmetry == Matrix_Symmetry::NONE || OutSymmetry == Matrix_Symmetry::STRUCTURALLY_SYMMETRIC) && OutShape != Matrix_Shape::FULL) {
        eslog::error("Submatrix: invalid combination of parameters - when unsymmetric, matrices cannot be upper nor lower\n");
        return;
    }

    constexpr bool is_T_real = std::is_same<T,float>::value || std::is_same<T,double>::value;
    constexpr bool is_T_complex = std::is_same<T,std::complex<float>>::value || std::is_same<T,std::complex<double>>::value;

    I out_rows = end_row - start_row;
    I out_cols = end_col - start_col;
    bool is_out_block_symmetric = (start_row == start_col && end_row == end_col);

    output.shape = OutShape;
    if(is_out_block_symmetric) {
        output.type = input.type;
    }
    else {
        if constexpr (is_T_real)    { output.type = Matrix_Type::REAL_NONSYMMETRIC; }
        if constexpr (is_T_complex) { output.type = Matrix_Type::COMPLEX_NONSYMMETRIC; }
    }

    if constexpr(DoTrans) { output.resize(out_cols, out_rows); }
    else { output.resize(out_rows, out_cols); }
    std::fill(output.vals, output.vals + output.nnz, T(0));

    constexpr bool do_copy                  = ((OutSymmetry == Matrix_Symmetry::NONE || OutSymmetry == Matrix_Symmetry::STRUCTURALLY_SYMMETRIC) && !DoTrans && !DoConj) || (OutSymmetry == Matrix_Symmetry::SYMMETRIC && !DoConj && (OutShape == Matrix_Shape::FULL || InShape == OutShape)) || (OutSymmetry == Matrix_Symmetry::HERMITIAN && DoTrans == DoConj && (OutShape == Matrix_Shape::FULL || InShape == OutShape));
    constexpr bool do_transpose_copy        = ((OutSymmetry == Matrix_Symmetry::NONE || OutSymmetry == Matrix_Symmetry::STRUCTURALLY_SYMMETRIC) &&  DoTrans && !DoConj) || (OutSymmetry == Matrix_Symmetry::SYMMETRIC && !DoConj &&                                    InShape != OutShape)  || (OutSymmetry == Matrix_Symmetry::HERMITIAN && DoTrans != DoConj &&                                    InShape != OutShape);
    constexpr bool do_copy_conjugated       = ((OutSymmetry == Matrix_Symmetry::NONE || OutSymmetry == Matrix_Symmetry::STRUCTURALLY_SYMMETRIC) && !DoTrans &&  DoConj) || (OutSymmetry == Matrix_Symmetry::SYMMETRIC &&  DoConj && (OutShape == Matrix_Shape::FULL || InShape == OutShape)) || (OutSymmetry == Matrix_Symmetry::HERMITIAN && DoTrans != DoConj && (OutShape == Matrix_Shape::FULL || InShape == OutShape));
    constexpr bool do_transposed_conjugated = ((OutSymmetry == Matrix_Symmetry::NONE || OutSymmetry == Matrix_Symmetry::STRUCTURALLY_SYMMETRIC) &&  DoTrans &&  DoConj) || (OutSymmetry == Matrix_Symmetry::SYMMETRIC &&  DoConj &&                                    InShape != OutShape)  || (OutSymmetry == Matrix_Symmetry::HERMITIAN && DoTrans == DoConj &&                                    InShape != OutShape);

    for(I out_r = 0; out_r < out_rows; ++out_r)
    {
        I in_r = out_r + start_row;
        I in_row_start_idx = input.rows[in_r];
        I in_row_end_idx = input.rows[in_r+1];

        I i = in_row_start_idx;
        while(i < in_row_end_idx && input.cols[i] < start_col) { i++; }
        for(; i < in_row_end_idx; i++) {
            I c = input.cols[i];
            if(c >= end_col) {
                break;
            }
            I out_c = c - start_col;

            T val;
            T conjval;
            T *out_val;
            T *out_val_trans;

            val = input.vals[i];
            if constexpr(do_copy_conjugated || do_transposed_conjugated) {
                if constexpr(is_T_real)    { conjval = val; }
                if constexpr(is_T_complex) { conjval = std::conj(val); }
            }
            if constexpr(do_copy || do_copy_conjugated) { out_val = output.vals + denseMatCalcIndexRowMajor<OutShape>(out_r, out_c, output.ncols); }
            if constexpr(do_transpose_copy || do_transposed_conjugated) { out_val_trans = output.vals + denseMatCalcIndexRowMajor<OutShape>(out_c, out_r, output.ncols); }

            if constexpr(do_copy) { *out_val = val; }
            if constexpr(do_transpose_copy) { *out_val_trans = val; }
            if constexpr(do_copy_conjugated) { *out_val = conjval; }
            if constexpr(do_transposed_conjugated) { *out_val_trans = conjval; }

            // since I use = instead of +=, duplicated operations on diagonal entries should not matter
        }
    }
}

template <typename T, typename I, bool DoTrans, bool DoConj, Matrix_Symmetry OutSymmetry, Matrix_Shape InShape>
static void _submatrix(const Matrix_CSR<T, I> &input, Matrix_Dense<T, I> &output, I start_row, I end_row, I start_col, I end_col, Matrix_Shape out_shape)
{
    switch(out_shape) {
    case Matrix_Shape::LOWER:
        _submatrix<T, I, DoTrans, DoConj, OutSymmetry, InShape, Matrix_Shape::LOWER>(input, output, start_row, end_row, start_col, end_col);
        break;
    case Matrix_Shape::FULL:
        _submatrix<T, I, DoTrans, DoConj, OutSymmetry, InShape, Matrix_Shape::FULL >(input, output, start_row, end_row, start_col, end_col);
        break;
    case Matrix_Shape::UPPER:
        _submatrix<T, I, DoTrans, DoConj, OutSymmetry, InShape, Matrix_Shape::UPPER>(input, output, start_row, end_row, start_col, end_col);
        break;
    }
}

template <typename T, typename I, bool DoTrans, bool DoConj, Matrix_Symmetry OutSymmetry>
static void _submatrix(const Matrix_CSR<T, I> &input, Matrix_Dense<T, I> &output, I start_row, I end_row, I start_col, I end_col, Matrix_Shape out_shape, Matrix_Shape in_shape)
{
    switch(in_shape) {
    case Matrix_Shape::LOWER:
        _submatrix<T, I, DoTrans, DoConj, OutSymmetry, Matrix_Shape::LOWER>(input, output, start_row, end_row, start_col, end_col, out_shape);
        break;
    case Matrix_Shape::FULL:
        _submatrix<T, I, DoTrans, DoConj, OutSymmetry, Matrix_Shape::FULL >(input, output, start_row, end_row, start_col, end_col, out_shape);
        break;
    case Matrix_Shape::UPPER:
        _submatrix<T, I, DoTrans, DoConj, OutSymmetry, Matrix_Shape::UPPER>(input, output, start_row, end_row, start_col, end_col, out_shape);
        break;
    }
}

template <typename T, typename I, bool DoTrans, bool DoConj>
static void _submatrix(const Matrix_CSR<T, I> &input, Matrix_Dense<T, I> &output, I start_row, I end_row, I start_col, I end_col, Matrix_Shape out_shape, Matrix_Shape in_shape, Matrix_Symmetry out_symmetry)
{
    switch(out_symmetry) {
    case Matrix_Symmetry::NONE:
        _submatrix<T, I, DoTrans, DoConj, Matrix_Symmetry::NONE                     >(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape);
        break;
    case Matrix_Symmetry::STRUCTURALLY_SYMMETRIC:
        _submatrix<T, I, DoTrans, DoConj, Matrix_Symmetry::STRUCTURALLY_SYMMETRIC>(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape);
        break;
    case Matrix_Symmetry::SYMMETRIC:
        _submatrix<T, I, DoTrans, DoConj, Matrix_Symmetry::SYMMETRIC             >(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape);
        break;
    case Matrix_Symmetry::HERMITIAN:
        _submatrix<T, I, DoTrans, DoConj, Matrix_Symmetry::HERMITIAN             >(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape);
        break;
    }
}

template <typename T, typename I, bool DoTrans>
static void _submatrix(const Matrix_CSR<T, I> &input, Matrix_Dense<T, I> &output, I start_row, I end_row, I start_col, I end_col, Matrix_Shape out_shape, Matrix_Shape in_shape, Matrix_Symmetry out_symmetry, bool conj)
{
    if(conj) {
        _submatrix<T, I, DoTrans, true >(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape, out_symmetry);
    }
    else {
        _submatrix<T, I, DoTrans, false>(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape, out_symmetry);
    }
}

template <typename T, typename I>
static void _submatrix(const Matrix_CSR<T, I> &input, Matrix_Dense<T, I> &output, I start_row, I end_row, I start_col, I end_col, Matrix_Shape out_shape, Matrix_Shape in_shape, Matrix_Symmetry out_symmetry, bool conj, bool trans)
{
    if(trans) {
        _submatrix<T, I, true >(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape, out_symmetry, conj);
    }
    else {
        _submatrix<T, I, false>(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape, out_symmetry, conj);
    }
}


template <typename T, typename I, bool DoConj>
static void _submatrix(const Matrix_CSR<T, I> &input, Matrix_CSR<T, I> &output, I start_row, I end_row, I start_col, I end_col)
{
    // start is inclusive, end is exclusive; all must be valid ranges in the context of this matrix
    // make sure all the desired values are present in the matrix (are not missing due to symmetry)

    constexpr bool is_T_real = std::is_same<T,float>::value || std::is_same<T,double>::value;
    constexpr bool is_T_complex = std::is_same<T,std::complex<float>>::value || std::is_same<T,std::complex<double>>::value;

    I out_rows = end_row - start_row;
    I out_cols = end_col - start_col;

    bool is_out_block_symmetric = (start_row == start_col && end_row == end_col);
    // bool is_output_symmetric = (is_out_block_symmetric && getSymmetry(input.type) == Matrix_Symmetry::SYMMETRIC);
    // bool is_output_hermitian = (is_out_block_symmetric && getSymmetry(input.type) == Matrix_Symmetry::HERMITIAN);

    if(is_out_block_symmetric) {
        output.type = input.type;
        output.shape = input.shape;
    }
    else {
        if constexpr (is_T_real)    { output.type = Matrix_Type::REAL_NONSYMMETRIC; }
        if constexpr (is_T_complex) { output.type = Matrix_Type::COMPLEX_NONSYMMETRIC; }
        output.shape = Matrix_Shape::FULL;
    }

    std::vector<I> colidx_starts(out_rows);
    std::vector<I> colidx_ends(out_rows);

    for(I out_r = 0; out_r < out_rows; out_r++) {
        I r = out_r + start_row;
        I row_start_idx = input.rows[r];
        I row_end_idx = input.rows[r+1];
        I i = row_start_idx;
        while(i < row_end_idx && input.cols[i] < start_col) { i++; }
        colidx_starts[out_r] = i;
        while(i < row_end_idx && input.cols[i] < end_col) { i++; }
        colidx_ends[out_r] = i;
    }

    I out_nnz = 0;
    for(I out_r = 0; out_r < out_rows; out_r++) {
        out_nnz += colidx_ends[out_r] - colidx_starts[out_r];
    }

    output.resize(out_rows, out_cols, out_nnz);

    I curr_idx = 0;
    for(I out_r = 0; out_r < out_rows; out_r++) {
        output.rows[out_r] = curr_idx;
        I colidx_start = colidx_starts[out_r];
        I colidx_end = colidx_ends[out_r];
        for(I i = colidx_start; i < colidx_end; i++)
        {
            output.cols[curr_idx] = input.cols[i] - start_col;
            if constexpr(is_T_complex && DoConj) { output.vals[curr_idx] = std::conj(input.vals[i]); }
            else { output.vals[curr_idx] = input.vals[i]; }
            curr_idx++;
        }
    }
    output.rows[out_rows] = curr_idx;
}


template <template <typename, typename, typename> class Matrix, typename T, typename I>
inline void SpBLAS<Matrix, T, I>::submatrix(Matrix_Dense<T, I> &output, I start_row, I end_row, I start_col, I end_col, bool trans, bool conj, bool output_force_full)
{
    submatrix(*matrix, output, start_row, end_row, start_col, end_col, trans, conj, output_force_full);
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
inline void SpBLAS<Matrix, T, I>::submatrix(Matrix_CSR<T, I> &output, I start_row, I end_row, I start_col, I end_col, bool trans, bool conj, bool output_force_full)
{
    submatrix(*matrix, output, start_row, end_row, start_col, end_col, trans, conj, output_force_full);
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
inline void SpBLAS<Matrix, T, I>::submatrix(const Matrix_CSR<T, I> &input, Matrix_Dense<T, I> &output, I start_row, I end_row, I start_col, I end_col, bool trans, bool conj, bool output_force_full)
{
    // I out_rows = end_row - start_row;
    // I out_cols = end_col - start_col;
    bool is_out_block_symmetric = (start_row == start_col && end_row == end_col);

    Matrix_Symmetry out_symmetry;
    if(is_out_block_symmetric) { out_symmetry = getSymmetry(input.type); }
    else { out_symmetry = Matrix_Symmetry::NONE; }

    Matrix_Shape in_shape = input.shape;

    Matrix_Shape out_shape;
    if(output_force_full || out_symmetry == Matrix_Symmetry::NONE || out_symmetry == Matrix_Symmetry::STRUCTURALLY_SYMMETRIC) { out_shape = Matrix_Shape::FULL; }
    else { out_shape = in_shape; }

    _submatrix<T, I>(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape, out_symmetry, conj, trans);
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
inline void SpBLAS<Matrix, T, I>::submatrix(const Matrix_CSR<T, I> &input, Matrix_CSR<T, I> &output, I start_row, I end_row, I start_col, I end_col, bool trans, bool conj, bool output_force_full)
{
    if(trans) {
        eslog::error("Submatrix CSR->CSR: transposition is not supported.\n");
    }
    if(output_force_full) {
        eslog::error("Submatrix CSR->CSR: forcing full output matrices is not supported.\n");
    }

    if(conj) {
        _submatrix<T, I, true >(input, output, start_row, end_row, start_col, end_col);
    }
    else {
        _submatrix<T, I, false>(input, output, start_row, end_row, start_col, end_col);
    }
}

}
