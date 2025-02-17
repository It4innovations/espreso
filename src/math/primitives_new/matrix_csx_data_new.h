
#ifndef SRC_MATH_PRIMITIVES_NEW_MATRIX_CSX_DATA_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MATRIX_CSX_DATA_NEW_H_

#include "math/primitives_new/matrix_csx_view_new.h"



template<typename T, typename I>
struct MatrixCsxData_new : public MatrixCsxView_new
{
    using MatrixCsxView_new::ator;
    using MatrixCsxView_new::ptrs;
    using MatrixCsxView_new::idxs;
    using MatrixCsxView_new::vals;
    using MatrixCsxView_new::nnz;
    using MatrixCsxView_new::order;
    using MatrixBase_new::nrows;
    using MatrixBase_new::ncols;
    using MatrixBase_new::diag;
    using MatrixBase_new::uplo;

    MatrixCsxData_new(const MatrixCsxData_new &) = delete;
    MatrixCsxData_new(MatrixCsxData_new && other)
    {
        std::swap(static_cast<MatrixCsxView_new&>(*this), static_cast<MatrixCsxView_new>(other));
    }
    MatrixCsxData_new & operator=(const MatrixCsxData_new &) = delete;
    MatrixCsxData_new & operator=(MatrixCsxData_new && other)
    {
        if(this == &other) return;
        std::swap(static_cast<MatrixCsxView_new&>(*this), static_cast<MatrixCsxView_new>(other));
        other.free();
    }

    MatrixCsxData_new()
    {
    }
    MatrixCsxData_new(size_t nrows_, size_t ncols_, size_t nnz_, char order_, Allocator_new * ator_ = AllocatorCPU_new::get_singleton())
    {
        realloc(nrows_, ncols_, nnz_, order_, ator_);
    }
    virtual ~MatrixCsxData_new()
    {
        free();
    }
    void resize(size_t nrows_, size_t ncols_, size_t nnz_, char order_, Allocator_new * ator_ = AllocatorCPU_new::get_singleton())
    {
        free();

        nrows = nrows_;
        ncols = ncols_;
        nnz = nnz_;
        order = order_;
        ator = ator_;

        alloc();
    }
    void alloc()
    {
        if(nnz > 0) {
            ptrs = ator->alloc<I>(get_primary_size() + 1);
            idxs = ator->alloc<I>(nnz);
            vals = ator->alloc<T>(nnz);
        }
    }
    void free()
    {
        if(ator != nullptr) {
            ator->free(ptrs);
            ator->free(idxs);
            ator->free(vals);
        }
    }
};



#endif /* SRC_MATH_PRIMITIVES_NEW_MATRIX_CSX_DATA_NEW_H_ */
