
#ifndef SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_DATA_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_DATA_NEW_H_

#include "math/primitives_new/matrix_dense_view_new.h"



namespace espreso {



template<typename T>
class MatrixDenseData_new : public MatrixDenseView_new<T>
{
public:
    using MatrixDenseView_new<T>::ator;
    using MatrixDenseView_new<T>::vals;
    using MatrixDenseView_new<T>::ld;
    using MatrixDenseView_new<T>::order;
    using MatrixDenseView_new<T>::was_set;
    using MatrixBase_new::nrows;
    using MatrixBase_new::ncols;
    using MatrixBase_new::prop;
public:
    MatrixDenseData_new() = default;
    MatrixDenseData_new(const MatrixDenseData_new &) = delete;
    MatrixDenseData_new(MatrixDenseData_new && other)
    {
        std::swap(*static_cast<MatrixDenseView_new<T>*>(this), *static_cast<MatrixDenseView_new<T>*>(&other));
    }
    MatrixDenseData_new & operator=(const MatrixDenseData_new &) = delete;
    MatrixDenseData_new & operator=(MatrixDenseData_new && other)
    {
        if(this != &other) {
            std::swap(*static_cast<MatrixDenseView_new<T>*>(this), *static_cast<MatrixDenseView_new<T>*>(&other));
            other.free();
        }
        return *this;
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
        size_t align = ator->get_align();
        if(align % sizeof(T) != 0) eslog::error("incompatible allocator alignment and sizeof(T)\n");
        align /= sizeof(T);
        ld = ((this->get_size_secdary() - 1) / align + 1) * align;
        was_set = true;
    }
    void alloc()
    {
        if(!was_set) eslog::error("matrix has not been set\n");
        if(ator == nullptr) eslog::error("matrix data has not been set\n");
        if(vals != nullptr) eslog::error("matrix is already allocated\n");
        if(nrows > 0 && ncols > 0) {
            vals = ator->template alloc<T>(this->get_size_primary() * ld);
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
public:
    size_t get_memory_impact() const
    {
        size_t mem_secdary = this->get_size_secdary() * sizeof(T);
        mem_secdary = utils::round_up(mem_secdary, ator->get_align());
        size_t mem_total = mem_secdary * this->get_size_primary();
        return mem_total;
    }
};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_MATRIX_DENSE_DATA_NEW_H_ */
