
#ifndef SRC_MATH_PRIMITIVES_NEW_MATRIX_CSX_VIEW_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MATRIX_CSX_VIEW_NEW_H_

#include "math/primitives_new/matrix_base_new.h"
#include "math/primitives_new/allocator_new.h"



template<typename T, typename I>
struct MatrixCsxView_new : public MatrixBase_new
{
    Allocator_new * ator;
    I * ptrs = nullptr;
    I * idxs = nullptr;
    T * vals = nullptr;
    size_t nnz = 0;
    char order = '_';
    using MatrixBase_new::nrows;
    using MatrixBase_new::ncols;
    using MatrixBase_new::diag;
    using MatrixBase_new::uplo;

    MatrixCsxView_new(const MatrixCsxView_new &) = default;
    MatrixCsxView_new(MatrixCsxView_new &&) = default;
    MatrixCsxView_new & operator=(const MatrixCsxView_new &) = default;
    MatrixCsxView_new & operator=(MatrixCsxView_new &&) = default;

protected:
    MatrixCsxView_new()
    {
    }
public:
    virtual ~MatrixCsxView_new()
    {
    }

    size_t get_primary_size()
    {
        if(order == 'R') return nrows;
        if(order == 'C') return ncols;
        return 0;
    }

    MatrixCsxView_new<T> make_transposed_reordered() const
    {
        MatrixCsxView_new M = *this;
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
    static MatrixCsxView_new<T> from_old(MatrixCSR<I,A> & M_old)
    {
        MatrixCsxView_new<T> M_new;
        M_new.ator = &AllocatorDummy_new::get_singleton(M_old::is_data_host_accessible, M_old::is_data_device_accessible);
        M_new.ptrs = M_old.rows;
        M_new.idxs = M_old.cols;
        M_new.vals = M_old.vals;
        M_new.nnz = M_old.nnz;
        M_new.order = 'R';
        M_new.nrows = M_old.nrows;
        M_new.ncols = M_old.ncols;
        M_new.diag = '_';
        M_new.uplo = '_';
        if(M_old.shape == Matrix_Shape::LOWER) M_new.uplo = 'L';
        if(M_old.shape == Matrix_Shape::UPPER) M_new.uplo = 'U';
        if(M_old.shape == Matrix_Shape::FULL) M_new.uplo = 'F';
        return M_new;
    }
    template<typename I, typename A>
    static MatrixCSR<I,A> to_old(MatrixCsxView_new<T> & M_new)
    {
        if(M_new.ator.is_on_cpu() != A::is_data_host_accessible || M_new.ator.is_on_gpu() != A::is_data_device_accessible) eslog::error("allocators not compatible\n");
        MatrixCSR<T,I,A> M_old;
        M_old.nrows = M_new.nrows;
        M_old.ncols = M_new.ncols;
        M_old.nnz = M_new.nnz;
        M_old.rows = M_new.ptrs;
        M_old.cols = M_new.idxs;
        M_old.vals = M_new.vals;
        if(M_new.uplo == 'U') M_old.shape == Matrix_Shape::UPPER;
        if(M_new.uplo == 'L') M_old.shape == Matrix_Shape::LOWER;
        if(M_new.uplo == 'F') M_old.shape == Matrix_Shape::FULL;
        return M_old;
    }
};



#endif /* SRC_MATH_PRIMITIVES_NEW_MATRIX_CSX_VIEW_NEW_H_ */
