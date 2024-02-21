
#ifndef SRC_MATH_MATH_H_
#define SRC_MATH_MATH_H_

#include "esinfo/eslog.h"
#include "primitives/vector_dense.h"
#include "primitives/vector_sparse.h"
#include "primitives/matrix_dense.h"
#include "primitives/matrix_csr.h"
#include "primitives/matrix_ijv.h"

#include "wrappers/math.blas.h"
#include "wrappers/math.spblas.h"
#include "wrappers/math.lapack.h"
#include "wrappers/math.spsolver.h"

#include "basis/utilities/utils.h"

#include <complex>
#include <vector>
#include <algorithm>

namespace espreso {
namespace math { // interface to wrappers

    // x = value
    template <typename T, typename I>
    void set(const I size, T *x, const int incX, const T &value)
    {
        for (esint i = 0; i < size; i += incX) {
            x[i] = value;
        }
    }

} // math (interface to wrappers)
} // espreso

namespace espreso {
namespace math {

    template <typename T, typename I> void copy(Vector_Dense<T, I>  &x, const Vector_Dense<T, I>  &y) { blas::copy(x.size           , x.vals, 1, y.vals, 1); }
    template <typename T, typename I> void copy(Vector_Sparse<T, I> &x, const Vector_Sparse<T, I> &y) { blas::copy(x.nnz            , x.vals, 1, y.vals, 1); }
    template <typename T, typename I> void copy(Matrix_Dense<T, I>  &x, const Matrix_Dense<T, I>  &y) { blas::copy(x.nrows * x.ncols, x.vals, 1, y.vals, 1); }
    template <typename T, typename I> void copy(Matrix_CSR<T, I>    &x, const Matrix_CSR<T, I>    &y) { blas::copy(x.nnz            , x.vals, 1, y.vals, 1); }
    template <typename T, typename I> void copy(Matrix_IJV<T, I>    &x, const Matrix_IJV<T, I>    &y) { blas::copy(x.nnz            , x.vals, 1, y.vals, 1); }

    template <typename T, typename I> void scale(const T &alpha, Vector_Dense<T, I>  &x) { blas::scale(x.size           , alpha, x.vals, 1); }
    template <typename T, typename I> void scale(const T &alpha, Vector_Sparse<T, I> &x) { blas::scale(x.nnz            , alpha, x.vals, 1); }
    template <typename T, typename I> void scale(const T &alpha, Matrix_Dense<T, I>  &x) { blas::scale(x.nrows * x.ncols, alpha, x.vals, 1); }
    template <typename T, typename I> void scale(const T &alpha, Matrix_CSR<T, I>    &x) { blas::scale(x.nnz            , alpha, x.vals, 1); }
    template <typename T, typename I> void scale(const T &alpha, Matrix_IJV<T, I>    &x) { blas::scale(x.nnz            , alpha, x.vals, 1); }

    template <typename T, typename I> void add(Vector_Dense<T, I>  &x, const T &alpha, const Vector_Dense<T, I>  &y) { blas::add(x.size           , x.vals, 1, alpha, y.vals, 1); }
    template <typename T, typename I> void add(Vector_Sparse<T, I> &x, const T &alpha, const Vector_Sparse<T, I> &y) { blas::add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }
    template <typename T, typename I> void add(Matrix_Dense<T, I>  &x, const T &alpha, const Matrix_Dense<T, I>  &y) { blas::add(x.nrows * x.ncols, x.vals, 1, alpha, y.vals, 1); }
    template <typename T, typename I> void add(Matrix_CSR<T, I>    &x, const T &alpha, const Matrix_CSR<T, I>    &y) { blas::add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }
    template <typename T, typename I> void add(Matrix_IJV<T, I>    &x, const T &alpha, const Matrix_IJV<T, I>    &y) { blas::add(x.nnz            , x.vals, 1, alpha, y.vals, 1); }

    template <typename T, typename I> T dot(const Vector_Dense<T, I>  &x, const Vector_Dense<T, I>  &y) { return blas::dot(x.size, x.vals, 1, y.vals, 1); }
    template <typename T, typename I> T dot(const Vector_Sparse<T, I> &x, const Vector_Sparse<T, I> &y) { return blas::dot(x.nnz , x.vals, 1, y.vals, 1); }

    template <typename T, typename I> T norm(const Vector_Dense<T, I>  &x) { return blas::norm(x.size, x.vals, 1); }
    template <typename T, typename I> T norm(const Vector_Sparse<T, I> &x) { return blas::norm(x.nnz , x.vals, 1); }

    template <typename T, typename I> void set(Vector_Dense<T, I>  &x, const T &value) { set(x.size           , x.vals, 1, value); }
    template <typename T, typename I> void set(Vector_Sparse<T, I> &x, const T &value) { set(x.nnz            , x.vals, 1, value); }
    template <typename T, typename I> void set(Matrix_Dense<T, I>  &x, const T &value) { set(x.nrows * x.ncols, x.vals, 1, value); }
    template <typename T, typename I> void set(Matrix_CSR<T, I>    &x, const T &value) { set(x.nnz            , x.vals, 1, value); }
    template <typename T, typename I> void set(Matrix_IJV<T, I>    &x, const T &value) { set(x.nnz            , x.vals, 1, value); }

    template <typename T, typename I> void combine(Matrix_CSR<T, I> &C, const Matrix_CSR<T, I> &A, const Matrix_CSR<T, I> &B)
    {
        if (A.nrows != B.nrows || A.ncols != B.ncols) {
            eslog::error("invalid matrices sizes.\n");
        }
//        if(A.type == B.type && A.shape == B.shape) {
//            C.type = A.type;
//            C.shape = A.shape;
//        }
//        else {
//            if(utils::is_real<T>()) C.type = Matrix_Type::REAL_NONSYMMETRIC;
//            else C.type = Matrix_Type::COMPLEX_NONSYMMETRIC;
//            C.shape = Matrix_Shape::FULL;
//        }
        esint nnz = 0;
        for (esint r = 0; r < A.nrows; ++r) {
            esint *beginA = A.cols + A.rows[r    ] - Indexing::CSR;
            esint *endA   = A.cols + A.rows[r + 1] - Indexing::CSR;
            esint *beginB = B.cols + B.rows[r    ] - Indexing::CSR;
            esint *endB   = B.cols + B.rows[r + 1] - Indexing::CSR;
            while (true) {
                if (beginA == endA) { nnz += endB - beginB; break; }
                if (beginB == endB) { nnz += endA - beginA; break; }
                if (*beginA == *beginB) {
                    ++beginA; ++beginB;
                } else {
                    if (*beginA < *beginB) {
                        ++beginA;
                    } else {
                        ++beginB;
                    }
                }
                ++nnz;
            }
        }
        C.resize(A.nrows, A.ncols, nnz);
        esint *r = C.rows, *c = C.cols;
        for (esint i = 0; i < A.nrows; ++i, ++r) {
            *r = c - C.cols + Indexing::CSR;
            esint *bA = A.cols + A.rows[i    ] - Indexing::CSR;
            esint *eA = A.cols + A.rows[i + 1] - Indexing::CSR;
            esint *bB = B.cols + B.rows[i    ] - Indexing::CSR;
            esint *eB = B.cols + B.rows[i + 1] - Indexing::CSR;
            while (true) {
                if (bA == eA) { while (bB != eB) { *c++ = *bB++;} break; }
                if (bB == eB) { while (bA != eA) { *c++ = *bA++;} break; }
                if (*bA == *bB) {
                    *c++ = *bA++; bB++;
                } else {
                    if (*bA < *bB) {
                        *c++ = *bA++;
                    } else {
                        *c++ = *bB++;
                    }
                }
            }
        }
        *r = c - C.cols + Indexing::CSR;
    }

    template <typename T, typename I> void sumCombined(Matrix_CSR<T, I> &C, const T &alpha, const Matrix_CSR<T, I> &A, const Matrix_CSR<T, I> &B)
    {
        if (A.nrows != B.nrows || A.ncols != B.ncols) {
            eslog::error("invalid matrices sizes.\n");
        }
        for (esint r = 0; r < A.nrows; ++r) {
            esint *bA = A.cols + A.rows[r    ] - Indexing::CSR;
            esint *eA = A.cols + A.rows[r + 1] - Indexing::CSR;
            esint *bB = B.cols + B.rows[r    ] - Indexing::CSR;
            esint *eB = B.cols + B.rows[r + 1] - Indexing::CSR;
            esint *bC = C.cols + C.rows[r    ] - Indexing::CSR;
            esint *eC = C.cols + C.rows[r + 1] - Indexing::CSR;
            T      *a = A.vals + A.rows[r    ] - Indexing::CSR;
            T      *b = B.vals + B.rows[r    ] - Indexing::CSR;
            T      *c = C.vals + C.rows[r    ] - Indexing::CSR;
            while (bC != eC) {
                *c = 0;
                if (bA != eA && *bC == *bA) { *c += *a++; ++bA; }
                if (bB != eB && *bC == *bB) { *c += *b++; ++bB; }
                ++bC; ++c;
            }
        }
    }

    template <typename T, typename I> void orthonormalize(Matrix_Dense<T, I> &m)
    {
        for (esint r = 0; r < m.nrows; ++r) {
            for (esint rr = 0; rr < r; ++rr) {
                T scale = math::blas::dot(m.ncols, m.vals + rr * m.ncols, 1, m.vals + r * m.ncols, 1) / math::blas::dot(m.ncols, m.vals + rr * m.ncols, 1, m.vals + rr * m.ncols, 1);
                math::blas::add(m.ncols, m.vals + r * m.ncols, 1, -scale, m.vals + rr * m.ncols, 1);
            }
            math::blas::scale(m.ncols, T{1.} / math::blas::norm(m.ncols, m.vals + r * m.ncols, 1), m.vals + r * m.ncols, 1);
        }
    }

    template <class T> void store(const T &x, const char* file);

    template<typename T, typename I, typename Ao, typename Ai, typename Ap>
    void permuteColumns(Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, const Permutation<I,Ap> & perm)
    {
        static_assert(Ao::is_data_host_accessible, "permuteColumns: the allocator does not provide host-accessible memory");
        static_assert(Ai::is_data_host_accessible, "permuteColumns: the allocator does not provide host-accessible memory");
        static_assert(Ap::is_data_host_accessible, "permuteColumns: the allocator does not provide host-accessible memory");
        if(output.nrows != input.nrows || output.ncols != input.ncols || output.nnz != input.nnz) eslog::error("permuteColumns: output matrix has wrong dimensions\n");
        if(perm.size != input.ncols) eslog::error("permuteColumns: permutation has wrong size\n");

        struct colval{ I col; T val; colval(I c, T v){ col = c; val = v;} };

        std::copy_n(input.rows, input.nrows + 1, output.rows);

        for(I r = 0; r < input.nrows; r++)
        {
            I start = input.rows[r];
            I end = input.rows[r+1];
            std::vector<colval> colvals_out;
            colvals_out.reserve(end - start);
            for(I i = start; i < end; i++)
            {
                T val = input.vals[i];
                I col_in = input.cols[i];
                I col_out = perm.src_to_dst[col_in];
                colvals_out.emplace_back(col_out, val);
            }
            std::sort(colvals_out.begin(), colvals_out.end(), [&](const colval & l, const colval & r){ return l.col < r.col; });
            for(size_t j = 0; j < colvals_out.size(); j++) { output.cols[start + j] = colvals_out[j].col; output.vals[start + j] = colvals_out[j].val; }
        }
    }

    template<typename T>
    T my_conj(T val)
    {
        if constexpr(utils::is_real<T>()) return val;
        if constexpr(utils::is_complex<T>()) return std::conj(val);
    }

    template<typename T, typename I, typename Ao, typename Am, typename Ai>
    void conjTransposeMapSetup(Matrix_CSR<T,I,Ao> & output, Vector_Dense<I,I,Am> & map, const Matrix_CSR<T,I,Ai> & input)
    {
        static_assert(Ao::is_data_host_accessible, "conjTransposeMapSetup: the allocator does not provide host-accessible memory");
        static_assert(Am::is_data_host_accessible, "conjTransposeMapSetup: the allocator does not provide host-accessible memory");
        static_assert(Ai::is_data_host_accessible, "conjTransposeMapSetup: the allocator does not provide host-accessible memory");
        if(output.nrows != input.ncols || output.ncols != input.nrows || output.nnz != input.nnz) eslog::error("conjTransposeMapSetup: output matrix has wrong dimensions\n");
        if(input.nnz != map.size) eslog::error("conjTransposeMapSetup: input map has wrong dimensions\n");

        struct colvalidx{ T val; I col; I idx; colvalidx(I c, T v, I i){ col = c; val = v; idx = i;} };

        std::vector<std::vector<colvalidx>> out_rows(output.nrows);

        for(I r = 0; r < input.nrows; r++)
        {
            I start = input.rows[r];
            I end = input.rows[r+1];
            for(I i = start; i < end; i++)
            {
                I c = input.cols[i];
                T v = input.vals[i];
                out_rows[c].emplace_back(r, v, i);
            }
        }

        I curr_idx = 0;
        for(I row_out = 0; row_out < output.nrows; row_out++)
        {
            output.rows[row_out] = curr_idx;
            std::vector<colvalidx> & data = out_rows[row_out];
            for(size_t i = 0; i < data.size(); i++)
            {
                output.cols[curr_idx] = data[i].col;
                output.vals[curr_idx] = my_conj(data[i].val);
                map.vals[curr_idx] = data[i].idx;
                curr_idx++;
            }
        }
        output.rows[output.nrows] = curr_idx;
    }

    template<typename T, typename I, typename Ao, typename Am, typename Ai>
    void conjTransposeMapUse(Matrix_CSR<T,I,Ao> & output, const Vector_Dense<I,I,Am> & map, const Matrix_CSR<T,I,Ai> & input)
    {
        static_assert(Ao::is_data_host_accessible, "conjTransposeMapUse: the allocator does not provide host-accessible memory");
        static_assert(Am::is_data_host_accessible, "conjTransposeMapUse: the allocator does not provide host-accessible memory");
        static_assert(Ai::is_data_host_accessible, "conjTransposeMapUse: the allocator does not provide host-accessible memory");
        if(output.nrows != input.ncols || output.ncols != input.nrows || output.nnz != input.nnz) eslog::error("conjTransposeMapUse: output matrix has wrong dimensions\n");
        if(input.nnz != map.size) eslog::error("conjTransposeMapUse: input map has wrong dimensions\n");

        for(I i = 0; i < map.size; i++)
        {
            output.vals[i] = my_conj(input.vals[map.vals[i]]);
        }
    }

} // math
} // espreso

#endif /* SRC_MATH_MATH_H_ */
