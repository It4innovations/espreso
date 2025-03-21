
#ifndef SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_VIEW_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_VIEW_NEW_H_

#include <cstring>

#include "math/primitives_new/matrix_base_new.h"
#include "math/primitives/matrix_dense.h"
#include "basis/utilities/utils.h"



namespace espreso {



template<typename T>
class MatrixDenseView_new : public MatrixBase_new
{
public: // the user promises not to modify these values (I don't want to implement getters everywhere)
    T * vals = nullptr;
    size_t ld = 0;
    char order = '_';
    bool was_set = false;
public:
    using MatrixBase_new::nrows;
    using MatrixBase_new::ncols;
    using MatrixBase_new::prop;
public:
    MatrixDenseView_new() = default;
    MatrixDenseView_new(const MatrixDenseView_new &) = default;
    MatrixDenseView_new(MatrixDenseView_new &&) = default;
    MatrixDenseView_new & operator=(const MatrixDenseView_new &) = default;
    MatrixDenseView_new & operator=(MatrixDenseView_new &&) = default;
    virtual ~MatrixDenseView_new() = default;
public:
    void set_view(size_t nrows_, size_t ncols_, size_t ld_, size_t order_, T * vals_)
    {
        if(was_set) eslog::error("can only set yet-uninitialized matrix view\n");
        nrows = nrows_;
        ncols = ncols_;
        ld = ld_;
        order = order_;
        vals = vals_;
        was_set = true;
    }
public:
    size_t get_size_primary() const
    {
        if(order == 'R') return nrows;
        if(order == 'C') return ncols;
        return 0;
    }
    size_t get_size_secdary() const
    {
        if(order == 'R') return ncols;
        if(order == 'C') return nrows;
        return 0;
    }
    size_t get_stride_row() const
    {
        if(order == 'R') return ld;
        if(order == 'C') return 1;
        return 0;
    }
    size_t get_stride_col() const
    {
        if(order == 'R') return 1;
        if(order == 'C') return ld;
        return 0;
    }

    MatrixDenseView_new<T> get_submatrix_view(size_t row_start, size_t row_end, size_t col_start, size_t col_end) const
    {
        if(row_start > row_end || row_end > nrows || col_start > col_end || col_end > ncols) eslog::error("wrong submatrix\n");

        MatrixDenseView_new<T> M = *this;
        M.nrows = row_end - row_start;
        M.ncols = col_end - col_start;
        M.vals = M.vals + row_start * M.get_stride_row() + col_start * M.get_stride_col();
        return M;
    }
    MatrixDenseView_new<T> get_transposed_reordered_view() const
    {
        MatrixDenseView_new<T> M = *this;
        M.transpose_reorder_inplace();
        return M;
    }
    void transpose_reorder_inplace()
    {
        std::swap(nrows, ncols);
        order = change_order(order);
        prop.uplo = change_uplo(prop.uplo);
    }

    static bool are_interchangable(const MatrixDenseView_new & A, const MatrixDenseView_new & B)
    {
        return (A.nrows == B.nrows) && (A.ncols == B.ncols) && (A.order == B.order) && (A.prop.uplo == B.prop.uplo) && (A.prop.diag == B.prop.diag);
    }

    template<typename I, typename A>
    static MatrixDenseView_new<T> from_old(const Matrix_Dense<T,I,A> & M_old, char order = 'R')
    {
        MatrixDenseView_new<T> M_new;
        M_new.set_view(M_old.nrows, M_old.ncols, M_old.get_ld(), order, M_old.vals);
        if(M_old.shape == Matrix_Shape::LOWER) M_new.prop.uplo = 'L';
        if(M_old.shape == Matrix_Shape::UPPER) M_new.prop.uplo = 'U';
        if(M_old.shape == Matrix_Shape::FULL) M_new.prop.uplo = 'F';
        return M_new;
    }
    template<typename I, typename A>
    static Matrix_Dense<T,I,A> to_old(MatrixDenseView_new<T> & M_new)
    {
        Matrix_Dense<T,I,A> M_old;
        M_old.nrows = M_new.nrows;
        M_old.ncols = M_new.ncols;
        M_old.set_ld(M_new.ld);
        M_old.vals = M_new.vals;
        if(M_new.prop.uplo == 'U') M_old.shape = Matrix_Shape::UPPER;
        else if(M_new.prop.uplo == 'L') M_old.shape = Matrix_Shape::LOWER;
        else M_old.shape = Matrix_Shape::FULL;
        return M_old;
    }

    void print(const char * name = "")
    {
        if constexpr(utils::is_real<T>()) {
            eslog::info("Dense matrix %s, size %zux%zu, ld %zu, order '%c', uplo '%c', diag '%c'\n", name, nrows, ncols, ld, order, prop.uplo, prop.diag);
            for(size_t r = 0; r < nrows; r++) {
                for(size_t c = 0; c < ncols; c++) {
                    if constexpr(std::is_floating_point_v<T>) {
                        double v = (double)vals[r * get_stride_row() + c * get_stride_col()];
                        char str[100];
                        snprintf(str, sizeof(str), "%+11.3e", v);
                        if(strstr(str, "nan") != nullptr) eslog::info("   nan      ");
                        else if(strstr(str, "inf") != nullptr) eslog::info("  %cinf      ", v > 0 ? '+' : '-');
                        else if(v == 0) eslog::info("   0        ");
                        else eslog::info(" %+11.3e", v);
                    }
                    if constexpr(std::is_integral_v<T>) {
                        long long v = (long long)vals[r * get_stride_row() + c * get_stride_col()];
                        eslog::info(" %+11lld", v);
                    }
                }
                eslog::info("\n");
            }
            fflush(stdout);
        }
        if constexpr(utils::is_complex<T>()) {
            eslog::error("matrix print not yet supported for complex matrices\n");
        }
    }
};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_VIEW_NEW_H_ */
