
#ifndef SRC_MATH_PRIMITIVES_NEW_VECTOR_DENSE_VIEW_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_VECTOR_DENSE_VIEW_NEW_H_

#include "math/primitives_new/vector_base_new.h"
#include "math/primitives/vector_dense.h"



namespace espreso {



template<typename T>
struct VectorDenseView_new : public VectorBase_new
{
public: // the user promises not to modify these values (I don't want to implement getters everywhere)
    T * vals = nullptr;
    bool was_set = false;
public:
    using VectorBase_new::size;
public:
    VectorDenseView_new() = default;
    VectorDenseView_new(const VectorDenseView_new &) = default;
    VectorDenseView_new(VectorDenseView_new &&) = default;
    VectorDenseView_new & operator=(const VectorDenseView_new &) = default;
    VectorDenseView_new & operator=(VectorDenseView_new &&) = default;
    virtual ~VectorDenseView_new() = default;
public:
    void set_view(size_t size_, T * vals_)
    {
        if(was_set) eslog::error("can only set yet-uninitialized vector view\n");
        size = size_;
        vals = vals_;
        was_set = true;
    }
public:
    static bool are_interchangable(VectorDenseView_new & A, VectorDenseView_new & B)
    {
        return (A.size == B.size);
    }

    template<typename I, typename A>
    static VectorDenseView_new<T> from_old(Vector_Dense<T,I,A> & V_old)
    {
        VectorDenseView_new<T> V_new;
        V_new.set_view(V_old.size, V_old.vals);
        return V_new;
    }
    template<typename I, typename A>
    static Vector_Dense<T,I,A> to_old(VectorDenseView_new<T> & V_new)
    {
        if(V_new.ator.is_on_cpu() != A::is_data_host_accessible || V_new.ator.is_on_gpu() != A::is_data_device_accessible) eslog::error("allocators not compatible\n");
        Vector_Dense<T,I,A> V_old;
        V_old.size = V_new.size;
        V_old.vals = V_new.vals;
        return V_old;
    }
};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_VECTOR_DENSE_VIEW_NEW_H_ */
