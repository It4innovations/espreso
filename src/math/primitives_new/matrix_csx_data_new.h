
#ifndef SRC_MATH_PRIMITIVES_NEW_MATRIX_CSX_DATA_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MATRIX_CSX_DATA_NEW_H_

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/allocator_new.h"



namespace espreso {



template<typename T, typename I>
struct MatrixCsxData_new : public MatrixCsxView_new
{
public: // the user promises not to modify these values (I don't want to implement getters everywhere)
    Allocator_new * ator = nullptr;
public:
    using MatrixCsxView_new::ptrs;
    using MatrixCsxView_new::idxs;
    using MatrixCsxView_new::vals;
    using MatrixCsxView_new::nnz;
    using MatrixCsxView_new::order;
    using MatrixBase_new::nrows;
    using MatrixBase_new::ncols;
    using MatrixBase_new::prop;
public:
    MatrixCsxData_new() = default;
    MatrixCsxData_new(const MatrixCsxData_new &) = delete;
    MatrixCsxData_new(MatrixCsxData_new && other)
    {
        std::swap(static_cast<MatrixCsxView_new&>(*this), static_cast<MatrixCsxView_new>(other));
        std::swap(ator, other.ator);
    }
    MatrixCsxData_new & operator=(const MatrixCsxData_new &) = delete;
    MatrixCsxData_new & operator=(MatrixCsxData_new && other)
    {
        if(this == &other) return;
        std::swap(static_cast<MatrixCsxView_new&>(*this), static_cast<MatrixCsxView_new>(other));
        std::swap(ator, other.ator);
        other.free();
    }
    virtual ~MatrixCsxData_new()
    {
        free();
    }
public:
    void set(size_t nrows_, size_t ncols_, size_t nnz_, char order_, Allocator_new * ator_)
    {
        if(was_set) eslog::error("can only set yet-uninitialized matrix view\n");
        nrows = nrows_;
        ncols = ncols_;
        nnz = nnz_;
        order = order_;
        ator = ator_;
        was_set = true;
    }
    void alloc()
    {
        if(!was_set) eslog::error("matrix has not been set\n");
        if(ator == nullptr) eslog::error("matrix data has not been set\n");
        if(vals != nullptr) eslog::error("matrix is already allocated\n");
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
    void clear()
    {
        free();
        nrows = 0;
        ncols = 0;
        ator = nullptr;
        was_set = false;
    }
};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_MATRIX_CSX_DATA_NEW_H_ */
