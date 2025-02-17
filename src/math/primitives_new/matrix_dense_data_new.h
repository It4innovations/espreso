
#ifndef SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_NEW_H_

#include "math/primitives_new/matrix_base_view_new.h"



template<typename T>
struct MatrixDenseData_new : public MatrixDenseView_new
{
    // std::vector<MatrixDenseView_new> views; ???
    using MatrixDenseView_new::ator;
    using MatrixDenseView_new::vals;
    using MatrixDenseView_new::ld;
    using MatrixDenseView_new::order;
    using MatrixBase_new::nrows;
    using MatrixBase_new::ncols;
    using MatrixBase_new::diag;
    using MatrixBase_new::uplo;

    MatrixDenseData_new(const MatrixDenseData_new &) = delete;
    MatrixDenseData_new(MatrixDenseData_new && other)
    {
        std::swap(static_cast<MatrixDenseView_new&>(*this), static_cast<MatrixDenseView_new>(other));
    }
    MatrixDenseData_new & operator=(const MatrixDenseData_new &) = delete;
    MatrixDenseData_new & operator=(MatrixDenseData_new && other)
    {
        if(this == &other) return;
        std::swap(static_cast<MatrixDenseView_new&>(*this), static_cast<MatrixDenseView_new>(other));
        other.free();
    }

    MatrixDenseData_new()
    {
    }
    MatrixDenseData_new(size_t nrows_, size_t ncols_, char order_, Allocator_new * ator_ = AllocatorCPU_new::get_singleton())
    {
        realloc(nrows_, ncols_, order_, ator_);
    }
    virtual ~MatrixDenseData_new()
    {
        free();
    }
    void resize(size_t nrows_, size_t ncols_, char order_, Allocator_new * ator_ = AllocatorCPU_new::get_singleton())
    {
        free();

        nrows = nrows_;
        ncols = ncols_;
        order = order_;
        ator = ator_;

        alloc();
    }
    void alloc()
    {
        if(nrows > 0 && ncols > 0) {
            vals = ator->alloc_2d<T>(get_num_blocks(), get_block_size(), ld);
        }
    }
    void free()
    {
        if(ator != nullptr) {
            ator->free(vals);
        }
    }
};



#endif /* SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_NEW_H_ */
