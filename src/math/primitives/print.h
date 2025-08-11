
#ifndef SRC_MATH_PRIMITIVES_PRINT_H
#define SRC_MATH_PRIMITIVES_PRINT_H

#include <cstdio>
#include <cstring>

#include "math/primitives/vector_dense.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "basis/utilities/utils.h"



namespace espreso {
namespace math {


    template<typename T, typename I, typename A>
    void print_vector(const Vector_Dense<T,I,A> & vec, const char * name = "")
    {
        static_assert(A::is_data_host_accessible);

        if constexpr(utils::is_real<T>()) {
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
        if constexpr(utils::is_complex<T>()) {
            eslog::error("matrix print not supported for complex matrices\n");
        }
    }

    template<typename T, typename I, typename A>
    void print_matrix_dense(const Matrix_Dense<T,I,A> & matrix, const char * name = "")
    {
        static_assert(A::is_data_host_accessible);

        if constexpr(utils::is_real<T>()) {
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
        if constexpr(utils::is_complex<T>()) {
            eslog::error("matrix print not supported for complex matrices\n");
        }
    }

    template<typename T, typename I, typename A>
    void print_matrix_csr_as_dense(const Matrix_CSR<T,I,A> & matrix, const char * name = "")
    {
        static_assert(A::is_data_host_accessible);

        if constexpr(utils::is_real<T>()) {
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
        if constexpr(utils::is_complex<T>()) {
            eslog::error("matrix print not supported for complex matrices\n");
        }
    }

    template<typename T, typename I, typename A>
    void print_matrix_csr_arrays(const Matrix_CSR<T,I,A> & matrix, const char * name = "")
    {
        static_assert(A::is_data_host_accessible);

        if constexpr(utils::is_real<T>()) {
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
        if constexpr(utils::is_complex<T>()) {
            eslog::error("matrix print not supported for complex matrices\n");
        }
    }

    template<typename T, typename I, typename A>
    void print_matrix_csr_by_rows(const Matrix_CSR<T,I,A> & matrix, const char * name = "")
    {
        static_assert(A::is_data_host_accessible);

        if constexpr(utils::is_real<T>()) {
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
        if constexpr(utils::is_complex<T>()) {
            eslog::error("matrix print not supported for complex matrices\n");
        }
    }

}
}



#endif /* SRC_MATH_PRIMITIVES_PRINT_H */
