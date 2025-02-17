
#ifndef SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_VIEW_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_VIEW_NEW_H_

#include "math/primitives_new/matrix_base_new.h"
#include "math/primitives_new/allocator_new.h"




template<typename T>
struct MatrixDenseView_new : public MatrixBase_new
{
    Allocator_new * ator;
    T * vals = nullptr;
    size_t ld = 0;
    char order = '_';
    using MatrixBase_new::nrows;
    using MatrixBase_new::ncols;
    using MatrixBase_new::diag;
    using MatrixBase_new::uplo;

    MatrixDenseView_new(const MatrixDenseView_new &) = default;
    MatrixDenseView_new(MatrixDenseView_new &&) = default;
    MatrixDenseView_new & operator=(const MatrixDenseView_new &) = default;
    MatrixDenseView_new & operator=(MatrixDenseView_new &&) = default;

protected:
    MatrixDenseView_new()
    {
    }
public:
    virtual ~MatrixDenseView_new()
    {
    }

    size_t get_num_blocks() const
    {
        if(order == 'R') return nrows;
        if(order == 'C') return ncols;
        return 0;
    }
    size_t get_block_size() const
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

    MatrixDenseView_new<T> make_submatrix(size_t row_start, size_t row_end, size_t col_start, size_t col_end) const
    {
        if(row_start > row_end || row_end > nrows || col_start > col_end || col_end > ncols) eslog::error("wrong submatrix");

        MatrixDenseView_new<T> M = *this;
        M.nrows = row_end - row_start;
        M.ncols = col_end - col_start;
        M.vals = M.vals + row_start * M.get_stride_row() + col_start * M.get_stride_col();
        return M;
    }
    MatrixDenseView_new<T> make_transposed_reordered() const
    {
        MatrixDenseView_new<T> M = *this;
        M.transpose_reorder_inplace();
        return M;
    }
    void transpose_reorder_inplace()
    {
        std::swap(nrows, ncols);
        order = change_order(order);
        uplo = change_uplo(uplo);
    }

    template<typename I, typename A>
    static MatrixDenseView_new<T> from_old(MatrixDense<T,I,A> & M_old, char order = 'R')
    {
        MatrixDenseView_new<T> M_new;
        M_new.ator = &AllocatorDummy_new::get_singleton(M_old::is_data_host_accessible, M_old::is_data_device_accessible);
        M_new.vals = M_old.vals;
        M_new.ld = M_old.get_ld();
        M_new.order = order;
        M_new.nrows = M_old.nrows;
        M_new.ncols = M_old.ncols;
        M_new.diag = '_';
        M_new.shape = '_';
        if(M_old.shape == Matrix_Shape::LOWER) M_new.uplo = 'L';
        if(M_old.shape == Matrix_Shape::UPPER) M_new.uplo = 'U';
        if(M_old.shape == Matrix_Shape::FULL) M_new.uplo = 'F';
        return M_new;
    }
    template<typename I, typename A>
    static MatrixDense<T,I,A> to_old(MatrixDenseView_new<T> & M_new)
    {
        if(M_new.ator.is_on_cpu() != A::is_data_host_accessible || M_new.ator.is_on_gpu() != A::is_data_device_accessible) eslog::error("allocators not compatible\n");
        MatrixDense<T,I,A> M_old;
        M_old.nrows = M_new.nrows;
        M_old.ncols = M_new.ncols;
        M_old.set_ld(M_new.ld);
        M_old.vals = M_new.vals;
        if(M_new.uplo == 'U') M_old.shape == Matrix_Shape::UPPER;
        if(M_new.uplo == 'L') M_old.shape == Matrix_Shape::LOWER;
        if(M_new.uplo == 'F') M_old.shape == Matrix_Shape::FULL;
        return M_old;
    }
};



#endif /* SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_VIEW_NEW_H_ */
