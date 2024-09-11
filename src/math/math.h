
#ifndef SRC_MATH_MATH_H_
#define SRC_MATH_MATH_H_

#include "esinfo/eslog.h"
#include "primitives/vector_dense.h"
#include "primitives/vector_sparse.h"
#include "primitives/matrix_dense.h"
#include "primitives/matrix_csr.h"
#include "primitives/matrix_ijv.h"

#include "wrappers/math.blas.h"
#include "wrappers/math.solver.h"
#include "wrappers/math.spblas.h"
#include "wrappers/math.lapack.h"
#include "wrappers/math.spsolver.h"

#include "basis/utilities/utils.h"

#include <complex>
#include <vector>
#include <algorithm>
#include <cstring>

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

    template <typename T, typename I> void add(Vector_Dense<T, I>  &x, const T &alpha, const Vector_Dense<T, I>  &y) { blas::add(x.size, x.vals, 1, alpha, y.vals, 1); }
    template <typename T, typename I> void add(Vector_Sparse<T, I> &x, const T &alpha, const Vector_Sparse<T, I> &y) { blas::add(x.nnz , x.vals, 1, alpha, y.vals, 1); }
    template <typename T, typename I> void add(Matrix_Dense<T, I>  &x, const T &alpha, const Matrix_Dense<T, I>  &y) { blas::add(x.nnz , x.vals, 1, alpha, y.vals, 1); }
    template <typename T, typename I> void add(Matrix_CSR<T, I>    &x, const T &alpha, const Matrix_CSR<T, I>    &y) { blas::add(x.nnz , x.vals, 1, alpha, y.vals, 1); }
    template <typename T, typename I> void add(Matrix_IJV<T, I>    &x, const T &alpha, const Matrix_IJV<T, I>    &y) { blas::add(x.nnz , x.vals, 1, alpha, y.vals, 1); }

    template <typename T, typename I> T dot(const Vector_Dense<T, I>  &x, const Vector_Dense<T, I>  &y) { return blas::dot(x.size, x.vals, 1, y.vals, 1); }
    template <typename T, typename I> T dot(const Vector_Sparse<T, I> &x, const Vector_Sparse<T, I> &y) { return blas::dot(x.nnz , x.vals, 1, y.vals, 1); }

    template <typename T, typename I> T norm(const Vector_Dense<T, I>  &x) { return blas::norm(x.size, x.vals, 1); }
    template <typename T, typename I> T norm(const Vector_Sparse<T, I> &x) { return blas::norm(x.nnz , x.vals, 1); }

    template <typename T, typename I> void set(Vector_Dense<T, I>  &x, const T &value) { set(x.size           , x.vals, 1, value); }
    template <typename T, typename I> void set(Vector_Sparse<T, I> &x, const T &value) { set(x.nnz            , x.vals, 1, value); }
    template <typename T, typename I> void set(Matrix_Dense<T, I>  &x, const T &value) { set(x.nrows * x.ncols, x.vals, 1, value); }
    template <typename T, typename I> void set(Matrix_CSR<T, I>    &x, const T &value) { set(x.nnz            , x.vals, 1, value); }
    template <typename T, typename I> void set(Matrix_IJV<T, I>    &x, const T &value) { set(x.nnz            , x.vals, 1, value); }

    template <typename T, typename I> void copy(Matrix_Dense<T, I>  &x, const Matrix_CSR<T, I>  &y)
    {
        x.resize(y.nrows, y.ncols);
        math::set(x, T{0});
        for (I r = 0; r < y.nrows; ++r) {
            for (I ci = y.rows[r]; ci < y.rows[r + 1]; ++ci) {
                I c = y.cols[ci - y.rows[0]] - y.rows[0];
                x.vals[r * y.ncols + c] = y.vals[ci - y.rows[0]];
                if (y.shape != Matrix_Shape::FULL) {
                    x.vals[c * y.ncols + r] = y.vals[ci - y.rows[0]];
                }
            }
        }
    }

    template <typename T, typename I> void copy(Matrix_CSR<T, I> &x, const Matrix_Dense<T, I> &y)
    {
        if (x.shape == Matrix_Shape::UPPER) {
            for (I r = 0; r < x.nrows; ++r) {
                for (I ci = x.rows[r]; ci < x.rows[r + 1]; ++ci) {
                    I c = x.cols[ci - Indexing::CSR] - Indexing::CSR;
                    x.vals[ci - Indexing::CSR] = y.vals[r * y.ncols + c];
                }
            }
        } else {
            eslog::error("Implement non-UPPER copy.\n");
        }
    }

    template <typename T, typename I> void eye(Matrix_Dense<T, I>  &x, const T &value)
    {
        math::set(x, 0.);
        for (int i = 0; i < std::max(x.nrows, x.ncols); ++i) {
            x.vals[i * x.ncols + i] = value;
        }
    }

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
        I nnz = 0;
        for (I r = 0; r < A.nrows; ++r) {
            I *beginA = A.cols + A.rows[r    ] - Indexing::CSR;
            I *endA   = A.cols + A.rows[r + 1] - Indexing::CSR;
            I *beginB = B.cols + B.rows[r    ] - Indexing::CSR;
            I *endB   = B.cols + B.rows[r + 1] - Indexing::CSR;
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
        I *r = C.rows, *c = C.cols;
        for (I i = 0; i < A.nrows; ++i, ++r) {
            *r = c - C.cols + Indexing::CSR;
           I *bA = A.cols + A.rows[i    ] - Indexing::CSR;
           I *eA = A.cols + A.rows[i + 1] - Indexing::CSR;
           I *bB = B.cols + B.rows[i    ] - Indexing::CSR;
           I *eB = B.cols + B.rows[i + 1] - Indexing::CSR;
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
        for (I r = 0; r < A.nrows; ++r) {
            I *bA = A.cols + A.rows[r    ] - Indexing::CSR;
            I *eA = A.cols + A.rows[r + 1] - Indexing::CSR;
            I *bB = B.cols + B.rows[r    ] - Indexing::CSR;
            I *eB = B.cols + B.rows[r + 1] - Indexing::CSR;
            I *bC = C.cols + C.rows[r    ] - Indexing::CSR;
            I *eC = C.cols + C.rows[r + 1] - Indexing::CSR;
            T  *a = A.vals + A.rows[r    ] - Indexing::CSR;
            T  *b = B.vals + B.rows[r    ] - Indexing::CSR;
            T  *c = C.vals + C.rows[r    ] - Indexing::CSR;
            while (bC != eC) {
                *c = 0;
                if (bA != eA && *bC == *bA) { *c += *a++; ++bA; }
                if (bB != eB && *bC == *bB) { *c += *b++; ++bB; }
                ++bC; ++c;
            }
        }
    }

    template <typename T, typename I> void orthonormalize(Matrix_Dense<T, I> &A);
    template <typename T, typename I> void permute(Matrix_CSR<T, I> &A, const std::vector<I> &perm);
    template <typename T, typename I> void getKernel(Matrix_CSR<T, I> &A, Matrix_Dense<T, I> &R, Matrix_CSR<T, I> &regMat, I maxDefect, I scSize);

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

        std::vector<I> out_nnz_in_rows(output.nrows, I{0});

        for(I i = 0; i < input.nnz; i++) out_nnz_in_rows[input.cols[i]]++;

        output.rows[0] = 0;
        for(I r = 0; r < output.nrows; r++) output.rows[r+1] = output.rows[r] + out_nnz_in_rows[r];

        std::vector<I> curr_out_row_idxs(output.nrows);
        std::copy_n(output.rows, output.nrows, curr_out_row_idxs.begin());

        for(I r_in = 0; r_in < input.nrows; r_in++)
        {
            I start = input.rows[r_in];
            I end = input.rows[r_in+1];
            for(I i = start; i < end; i++)
            {
                I c_in = input.cols[i];
                I out_index = curr_out_row_idxs[c_in];
                curr_out_row_idxs[c_in]++;
                output.cols[out_index] = r_in;
                output.vals[out_index] = my_conj(input.vals[i]);
                map.vals[out_index] = i;
            }
        }
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

    template<typename T, typename I, typename A>
    void print_vector(Vector_Dense<T,I,A> & vec, const char * name = "")
    {
        static_assert(A::is_data_host_accessible);

        eslog::info("Vector %s, size %lld\n", name, (long long)vec.size);
        for(I i = 0; i < vec.size; i++) {
            if constexpr(std::is_floating_point_v<T>) {
                double v = (double)vec.vals[i];
                char str[100];
                snprintf(str, sizeof(str), "%+11.3e", v);
                if(strstr(str, "nan") != nullptr) eslog::info("   nan      ");
                else if(strstr(str, "inf") != nullptr) eslog::info("  %cinf      ", v > 0 ? '+' : '-');
                else if(v == 0) eslog::info("   0        ");
                else eslog::info(" %+11.3e", v);
            }
            if constexpr(std::is_integral_v<T>) eslog::info(" %+11lld", (long long)vec.vals[i]);
        }
        eslog::info("\n");
        fflush(stdout);
    }

    template<typename T, typename I, typename A>
    void print_matrix_dense(Matrix_Dense<T,I,A> & matrix, const char * name = "")
    {
        static_assert(A::is_data_host_accessible);

        eslog::info("Dense matrix %s, size %lldx%lld, ld %lld\n", name, (long long)matrix.nrows, (long long)matrix.ncols, (long long)matrix.get_ld());
        for(I r = 0; r < matrix.nrows; r++) {
            for(I c = 0; c < matrix.ncols; c++) {
                if constexpr(std::is_floating_point_v<T>) {
                    double v = (double)matrix.vals[r * matrix.get_ld() + c];
                    char str[100];
                    snprintf(str, sizeof(str), "%+11.3e", v);
                    if(strstr(str, "nan") != nullptr) eslog::info("   nan      ");
                    else if(strstr(str, "inf") != nullptr) eslog::info("  %cinf      ", v > 0 ? '+' : '-');
                    else if(v == 0) eslog::info("   0        ");
                    else eslog::info(" %+11.3e", v);
                }
                if constexpr(std::is_integral_v<T>) eslog::info(" %+11lld", (long long)matrix.vals[r * matrix.get_ld() + c]);
            }
            eslog::info("\n");
        }
        fflush(stdout);
    }

    template<typename T, typename I, typename A>
    void print_matrix_csr_as_dense(Matrix_CSR<T,I,A> & matrix, const char * name = "")
    {
        static_assert(A::is_data_host_accessible);

        eslog::info("CSR matrix %s, size %lldx%lld, nnz %lld\n", name, (long long)matrix.nrows, (long long)matrix.ncols, (long long)matrix.nnz);
        for(I r = 0; r < matrix.nrows; r++) {
            I start = matrix.rows[r];
            I end = matrix.rows[r+1];
            I curr_col = 0;
            for(I i = start; i <= end; i++) {
                I col = matrix.ncols;
                if(i < end) col = matrix.cols[i];
                for(I c = curr_col; c < col; c++) {
                    eslog::info("    .       ");
                }
                if(i < end) {
                    double val = matrix.vals[i];
                    if(val == 0) eslog::info("    0       ");
                    else eslog::info(" %+11.3e", val);
                    curr_col = col + 1;
                }
            }
            eslog::info("\n");
        }
        fflush(stdout);
    }

    template<typename T, typename I, typename A>
    void print_matrix_csr_arrays(Matrix_CSR<T,I,A> & matrix, const char * name = "")
    {
        static_assert(A::is_data_host_accessible);

        eslog::info("CSR matrix %s, size %lldx%lld, nnz %lld\n", name, (long long)matrix.nrows, (long long)matrix.ncols, (long long)matrix.nnz);
        eslog::info("row ptrs: ");
        for(I r = 0; r <= matrix.nrows; r++) eslog::info("%lld ", (long long)matrix.rows[r]);
        eslog::info("\n");
        eslog::info("col idxs: ");
        for(I i = 0; i < matrix.nnz; i++) eslog::info("%lld ", (long long)matrix.cols[i]);
        eslog::info("\n");
        eslog::info("vals:     ");
        for(I i = 0; i < matrix.nnz; i++) eslog::info("%+.3e ", (double)matrix.vals[i]);
        eslog::info("\n");
        fflush(stdout);
    }

    template<typename T, typename I, typename A>
    void print_matrix_csr_by_rows(Matrix_CSR<T,I,A> & matrix, const char * name = "")
    {
        static_assert(A::is_data_host_accessible);

        eslog::info("CSR matrix %s, size %lldx%lld, nnz %lld\n", name, (long long)matrix.nrows, (long long)matrix.ncols, (long long)matrix.nnz);
        for(I r = 0; r <= matrix.nrows; r++)
        {
            I start = matrix.rows[r];
            I end = matrix.rows[r+1];
            eslog::info("row %lld, indexes %lld--%lld:\n", (long long)r, (long long)start, (long long)end);
            eslog::info("colidxs:");
            for(I i = start; i < end; i++) eslog::info("%12lld, ", (long long)matrix.cols[i]);
            eslog::info("\n");
            eslog::info("vals:   ");
            for(I i = start; i < end; i++) eslog::info("%+12.3e, ", (double)matrix.vals[i]);
            eslog::info("\n");
        }
        fflush(stdout);
    }

} // math
} // espreso

#endif /* SRC_MATH_MATH_H_ */
