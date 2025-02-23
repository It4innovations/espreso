
#ifndef SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_NEW_H_

#include "math/primitives_new/matrix_base_view_new.h"
#include "math/primitives_new/allocator_new.h"



template<typename T>
class MatrixDenseData_new : public MatrixDenseView_new<T>
{
public: // the user promises not to modify these values (I don't want to implement getters everywhere)
    Allocator_new * ator = nullptr;
public:
    using MatrixDenseView_new::vals;
    using MatrixDenseView_new::ld;
    using MatrixDenseView_new::order;
    using MatrixDenseView_new::was_set;
    using MatrixBase_new::nrows;
    using MatrixBase_new::ncols;
    using MatrixBase_new::prop;
public:
    MatrixDenseData_new() = default;
    MatrixDenseData_new(const MatrixDenseData_new &) = delete;
    MatrixDenseData_new(MatrixDenseData_new && other)
    {
        std::swap(static_cast<MatrixDenseView_new&>(*this), static_cast<MatrixDenseView_new&>(other));
        std::swap(ator, other.ator);
    }
    MatrixDenseData_new & operator=(const MatrixDenseData_new &) = delete;
    MatrixDenseData_new & operator=(MatrixDenseData_new && other)
    {
        if(this == &other) return;
        std::swap(static_cast<MatrixDenseView_new&>(*this), static_cast<MatrixDenseView_new&>(other));
        std::swap(ator, other.ator);
        other.free();
    }
    virtual ~MatrixDenseData_new()
    {
        free();
    }
public:
    void set(size_t nrows_, size_t ncols_, char order_, Allocator_new * ator_)
    {
        if(was_set) eslog::error("can only set yet-uninitialized matrix data\n");
        nrows = nrows_;
        ncols = ncols_;
        order = order_;
        ator = ator_;
        was_set = true;
    }
    void alloc()
    {
        if(!was_set) eslog::error("matrix has not been set\n");
        if(ator == nullptr) eslog::error("matrix data has not been set\n");
        if(vals != nullptr) eslog::error("matrix is already allocated\n");
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
    void clear()
    {
        free();
        nrows = 0;
        ncols = 0;
        ator = nullptr;
        was_set = false;
    }
};



#endif /* SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_NEW_H_ */
