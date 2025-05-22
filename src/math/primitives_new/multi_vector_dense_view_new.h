
#ifndef SRC_MATH_PRIMITIVES_NEW_MULTI_VECTOR_DENSE_VIEW_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MULTI_VECTOR_DENSE_VIEW_NEW_H_

#include "math/primitives_new/vector_base_new.h"
#include "math/primitives_new/allocator_new_base.h"
#include "basis/utilities/utils.h"

#include <cstring>



namespace espreso {



template<typename T, typename I>
struct MultiVectorDenseView_new : public VectorBase_new
{
public: // the user promises not to modify these values (I don't want to implement getters everywhere)
    Allocator_new * ator = nullptr;
    size_t num_vectors = 0;
    I * offsets = nullptr;
    T * vals = nullptr;
    bool was_set = false;
public:
    using VectorBase_new::size;
public:
    MultiVectorDenseView_new() = default;
    MultiVectorDenseView_new(const MultiVectorDenseView_new &) = default;
    MultiVectorDenseView_new(MultiVectorDenseView_new &&) = default;
    MultiVectorDenseView_new & operator=(const MultiVectorDenseView_new &) = default;
    MultiVectorDenseView_new & operator=(MultiVectorDenseView_new &&) = default;
    virtual ~MultiVectorDenseView_new() = default;
public:
    void set_view(size_t num_vectors_, size_t total_size_, I * offsets_, T * vals_, Allocator_new * ator_)
    {
        if(was_set) eslog::error("can only set yet-uninitialized multivector view\n");
        num_vectors = num_vectors_;
        size = total_size_;
        offsets = offsets_;
        vals = vals_;
        ator = ator_;
        was_set = true;
    }
    T & at(size_t idx_outer, size_t idx_inner)
    {
        // assume the user knows what they are doing and does not access invalid memory, e.g. in GPU
        return vals[offsets[idx_outer] + idx_inner];
    }
    // cannot have get_vector_view(i), because offsets might be in GPU memory which is inaccessible, and I need to offsets to get the size of i-th vector
};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_MULTI_VECTOR_DENSE_VIEW_NEW_H_ */
