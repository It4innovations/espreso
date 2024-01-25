
#include "math/wrappers/math.spblas.h"
#include "math/wrappers/math.spsolver.h"
#include "esinfo/eslog.h"

#include <algorithm>
#include <complex>
#include <type_traits>
#include <vector>

namespace espreso {

#ifndef HAVE_MKL
#ifndef HAVE_SUITESPARSE



template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::SpBLAS()
: matrix{}, _spblas{}
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::~SpBLAS()
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
SpBLAS<Matrix, T, I>::SpBLAS(const MatrixType &a)
: matrix{}, _spblas{}
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::insert(const MatrixType &a)
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::insertTransposed(const MatrixType &a)
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::extractUpper(MatrixType &a)
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::apply(Vector_Dense<T> &y, const T &alpha, const T &beta, const Vector_Dense<T> &x)
{
	eslog::error("calling of empty SpBLAS wrapper.\n");
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::transposeTo(SpBLAS<Matrix, T, I> &A)
{

}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::multiply(SpBLAS<Matrix, T, I> &A, SpBLAS<Matrix, T, I> &B)
{
	eslog::error("calling of empty SpBLAS wrapper.\n");
}

#endif
#endif

template<Matrix_Shape Shape>
static esint denseMatCalcIndexRowMajor(esint r, esint c, esint ncols)
{
	if constexpr(Shape == Matrix_Shape::LOWER) { return r * ncols + (r * (r - 1) / 2)     + c; }
	if constexpr(Shape == Matrix_Shape::FULL)  { return r * ncols                         + c; }
	if constexpr(Shape == Matrix_Shape::UPPER) { return r * ncols - (r * (r - 1) / 2) - r + c; }
	return 0;
}

template <typename T, bool DoTrans, bool DoConj, Matrix_Symmetry OutSymmetry, Matrix_Shape InShape, Matrix_Shape OutShape>
static void _submatrix(const Matrix_CSR<T> &input, Matrix_Dense<T> &output, esint start_row, esint end_row, esint start_col, esint end_col)
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

	esint out_rows = end_row - start_row;
	esint out_cols = end_col - start_col;
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

	for(esint out_r = 0; out_r < out_rows; ++out_r)
	{
		esint in_r = out_r + start_row;
		esint in_row_start_idx = input.rows[in_r];
		esint in_row_end_idx = input.rows[in_r+1];

		esint i = in_row_start_idx;
		while(i < in_row_end_idx && input.cols[i] < start_col) { i++; }
		for(; i < in_row_end_idx; i++) {
			esint c = input.cols[i];
			if(c >= end_col) {
				break;
			}
			esint out_c = c - start_col;

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

template <typename T, bool DoTrans, bool DoConj, Matrix_Symmetry OutSymmetry, Matrix_Shape InShape>
static void _submatrix(const Matrix_CSR<T> &input, Matrix_Dense<T> &output, esint start_row, esint end_row, esint start_col, esint end_col, Matrix_Shape out_shape)
{
	switch(out_shape) {
	case Matrix_Shape::LOWER:
		_submatrix<T, DoTrans, DoConj, OutSymmetry, InShape, Matrix_Shape::LOWER>(input, output, start_row, end_row, start_col, end_col);
		break;
	case Matrix_Shape::FULL:
		_submatrix<T, DoTrans, DoConj, OutSymmetry, InShape, Matrix_Shape::FULL >(input, output, start_row, end_row, start_col, end_col);
		break;
	case Matrix_Shape::UPPER:
		_submatrix<T, DoTrans, DoConj, OutSymmetry, InShape, Matrix_Shape::UPPER>(input, output, start_row, end_row, start_col, end_col);
		break;
	}
}

template <typename T, bool DoTrans, bool DoConj, Matrix_Symmetry OutSymmetry>
static void _submatrix(const Matrix_CSR<T> &input, Matrix_Dense<T> &output, esint start_row, esint end_row, esint start_col, esint end_col, Matrix_Shape out_shape, Matrix_Shape in_shape)
{
	switch(in_shape) {
	case Matrix_Shape::LOWER:
		_submatrix<T, DoTrans, DoConj, OutSymmetry, Matrix_Shape::LOWER>(input, output, start_row, end_row, start_col, end_col, out_shape);
		break;
	case Matrix_Shape::FULL:
		_submatrix<T, DoTrans, DoConj, OutSymmetry, Matrix_Shape::FULL >(input, output, start_row, end_row, start_col, end_col, out_shape);
		break;
	case Matrix_Shape::UPPER:
		_submatrix<T, DoTrans, DoConj, OutSymmetry, Matrix_Shape::UPPER>(input, output, start_row, end_row, start_col, end_col, out_shape);
		break;
	}
}

template <typename T, bool DoTrans, bool DoConj>
static void _submatrix(const Matrix_CSR<T> &input, Matrix_Dense<T> &output, esint start_row, esint end_row, esint start_col, esint end_col, Matrix_Shape out_shape, Matrix_Shape in_shape, Matrix_Symmetry out_symmetry)
{
	switch(out_symmetry) {
	case Matrix_Symmetry::NONE:
		_submatrix<T, DoTrans, DoConj, Matrix_Symmetry::NONE					 >(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape);
		break;
	case Matrix_Symmetry::STRUCTURALLY_SYMMETRIC:
		_submatrix<T, DoTrans, DoConj, Matrix_Symmetry::STRUCTURALLY_SYMMETRIC>(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape);
		break;
	case Matrix_Symmetry::SYMMETRIC:
		_submatrix<T, DoTrans, DoConj, Matrix_Symmetry::SYMMETRIC			 >(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape);
		break;
	case Matrix_Symmetry::HERMITIAN:
		_submatrix<T, DoTrans, DoConj, Matrix_Symmetry::HERMITIAN			 >(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape);
		break;
	}
}

template <typename T, bool DoTrans>
static void _submatrix(const Matrix_CSR<T> &input, Matrix_Dense<T> &output, esint start_row, esint end_row, esint start_col, esint end_col, Matrix_Shape out_shape, Matrix_Shape in_shape, Matrix_Symmetry out_symmetry, bool conj)
{
	if(conj) {
		_submatrix<T, DoTrans, true >(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape, out_symmetry);
	}
	else {
		_submatrix<T, DoTrans, false>(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape, out_symmetry);
	}
}

template <typename T>
static void _submatrix(const Matrix_CSR<T> &input, Matrix_Dense<T> &output, esint start_row, esint end_row, esint start_col, esint end_col, Matrix_Shape out_shape, Matrix_Shape in_shape, Matrix_Symmetry out_symmetry, bool conj, bool trans)
{
	if(trans) {
		_submatrix<T, true >(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape, out_symmetry, conj);
	}
	else {
		_submatrix<T, false>(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape, out_symmetry, conj);
	}
}


template <typename T, bool DoConj>
static void _submatrix(const Matrix_CSR<T> &input, Matrix_CSR<T> &output, esint start_row, esint end_row, esint start_col, esint end_col)
{
	// start is inclusive, end is exclusive; all must be valid ranges in the context of this matrix
	// make sure all the desired values are present in the matrix (are not missing due to symmetry)

	constexpr bool is_T_real = std::is_same<T,float>::value || std::is_same<T,double>::value;
	constexpr bool is_T_complex = std::is_same<T,std::complex<float>>::value || std::is_same<T,std::complex<double>>::value;

	esint out_rows = end_row - start_row;
	esint out_cols = end_col - start_col;

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

	std::vector<esint> colidx_starts(out_rows);
	std::vector<esint> colidx_ends(out_rows);

	for(esint out_r = 0; out_r < out_rows; out_r++) {
		esint r = out_r + start_row;
		esint row_start_idx = input.rows[r];
		esint row_end_idx = input.rows[r+1];
		esint i = row_start_idx;
		while(i < row_end_idx && input.cols[i] < start_col) { i++; }
		colidx_starts[out_r] = i;
		while(i < row_end_idx && input.cols[i] < end_col) { i++; }
		colidx_ends[out_r] = i;
	}

	esint out_nnz = 0;
	for(esint out_r = 0; out_r < out_rows; out_r++) {
		out_nnz += colidx_ends[out_r] - colidx_starts[out_r];
	}

	output.resize(out_rows, out_cols, out_nnz);

	esint curr_idx = 0;
	for(esint out_r = 0; out_r < out_rows; out_r++) {
		output.rows[out_r] = curr_idx;
		esint colidx_start = colidx_starts[out_r];
		esint colidx_end = colidx_ends[out_r];
		for(esint i = colidx_start; i < colidx_end; i++)
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
void SpBLAS<Matrix, T, I>::submatrix(Matrix_Dense<T, I> &output, esint start_row, esint end_row, esint start_col, esint end_col, bool trans, bool conj, bool output_force_full)
{
	submatrix(*matrix, output, start_row, end_row, start_col, end_col, trans, conj, output_force_full);
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::submatrix(Matrix_CSR<T, I> &output, esint start_row, esint end_row, esint start_col, esint end_col, bool trans, bool conj, bool output_force_full)
{
	submatrix(*matrix, output, start_row, end_row, start_col, end_col, trans, conj, output_force_full);
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::submatrix(const Matrix_CSR<T, I> &input, Matrix_Dense<T, I> &output, esint start_row, esint end_row, esint start_col, esint end_col, bool trans, bool conj, bool output_force_full)
{
	// esint out_rows = end_row - start_row;
	// esint out_cols = end_col - start_col;
	bool is_out_block_symmetric = (start_row == start_col && end_row == end_col);

	Matrix_Symmetry out_symmetry;
	if(is_out_block_symmetric) { out_symmetry = getSymmetry(input.type); }
	else { out_symmetry = Matrix_Symmetry::NONE; }

	Matrix_Shape in_shape = input.shape;

	Matrix_Shape out_shape;
	if(output_force_full || out_symmetry == Matrix_Symmetry::NONE || out_symmetry == Matrix_Symmetry::STRUCTURALLY_SYMMETRIC) { out_shape = Matrix_Shape::FULL; }
	else { out_shape = in_shape; }

	_submatrix<T>(input, output, start_row, end_row, start_col, end_col, out_shape, in_shape, out_symmetry, conj, trans);
}

template <template <typename, typename, typename> class Matrix, typename T, typename I>
void SpBLAS<Matrix, T, I>::submatrix(const Matrix_CSR<T, I> &input, Matrix_CSR<T, I> &output, esint start_row, esint end_row, esint start_col, esint end_col, bool trans, bool conj, bool output_force_full)
{
	if(trans) {
		eslog::error("Submatrix CSR->CSR: transposition is not supported.\n");
	}
	if(output_force_full) {
		eslog::error("Submatrix CSR->CSR: forcing full output matrices is not supported.\n");
	}

	if(conj) {
		_submatrix<T, true >(input, output, start_row, end_row, start_col, end_col);
	}
	else {
		_submatrix<T, false>(input, output, start_row, end_row, start_col, end_col);
	}
}

template struct SpBLAS<Matrix_CSR, float>;
template struct SpBLAS<Matrix_CSR, double>;
template struct SpBLAS<Matrix_CSR, std::complex<double>>;

}
