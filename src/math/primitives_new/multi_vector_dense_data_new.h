
#ifndef SRC_MATH_PRIMITIVES_NEW_MULTI_VECTOR_DENSE_DATA_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MULTI_VECTOR_DENSE_DATA_NEW_H_

#include "math/primitives_new/multi_vector_dense_view_new.h"

#include "basis/utilities/utils.h"



namespace espreso {



template<typename T, typename I>
struct MultiVectorDenseData_new : public MultiVectorDenseView_new<T,I>
{
public:
    using MultiVectorDenseView_new<T,I>::ator;
    using MultiVectorDenseView_new<T,I>::num_vectors;
    using MultiVectorDenseView_new<T,I>::offsets;
    using MultiVectorDenseView_new<T,I>::vals;
    using MultiVectorDenseView_new<T,I>::was_set;
    using VectorBase_new::size;
public:
    MultiVectorDenseData_new() = default;
    MultiVectorDenseData_new(const MultiVectorDenseData_new &) = delete;
    MultiVectorDenseData_new(MultiVectorDenseData_new && other)
    {
        std::swap(*static_cast<MultiVectorDenseView_new<T,I>*>(this), *static_cast<MultiVectorDenseView_new<T,I>*>(&other));
    }
    MultiVectorDenseData_new & operator=(const MultiVectorDenseData_new &) = delete;
    MultiVectorDenseData_new & operator=(MultiVectorDenseData_new && other)
    {
        if(this != &other) {
            std::swap(*static_cast<MultiVectorDenseView_new<T,I>*>(this), *static_cast<MultiVectorDenseView_new<T,I>*>(&other));
            other.free();
        }
        return *this;
    }
    virtual ~MultiVectorDenseData_new()
    {
        free();
    }
public:
    void set(size_t num_vectors_, size_t total_size_, Allocator_new * ator_)
    {
        if(was_set) eslog::error("can only set yet-uninitialized multivector data\n");
        num_vectors = num_vectors_;
        size = total_size_;
        ator = ator_;
        was_set = true;
    }
    void alloc()
    {
        if(!was_set) eslog::error("multivector has not been set\n");
        if(ator == nullptr) eslog::error("multivector has not been set\n");
        if(vals != nullptr) eslog::error("multivector already contains data\n");
        if(size > 0) {
            vals = ator->template alloc<T>(size);
        }
        offsets = ator->template alloc<I>(num_vectors + 1);
    }
    void free()
    {
        if(ator != nullptr) {
            ator->free(vals);
            ator->free(offsets);
        }
    }
    void clear()
    {
        free();
        num_vectors = 0;
        size = 0;
        ator = nullptr;
        was_set = false;
    }
public:
    static MultiVectorDenseData_new<T,I> convert_from(std::vector<std::vector<T>> & input, Allocator_new * ator)
    {
        if(!ator->is_data_accessible_cpu()) eslog::error("wrong allocator\n");

        size_t total_size = 0;
        for(size_t i = 0; i < input.size(); i++) total_size += input[i].size();

        if(total_size >= std::numeric_limits<I>::max()) eslog::error("insufficient integer range\n");

        MultiVectorDenseData_new<T,I> output;
        output.set(input.size(), total_size, ator);
        output.alloc();

        I curr_idx = 0;
        for(size_t i = 0; i < input.size(); i++) {
            output.offsets[i] = curr_idx;
            std::copy_n(input[i].data(), input[i].size(), output.vals + curr_idx);
            curr_idx += input[i].size();
        }
        output.offsets[input.size()] = curr_idx;

        return output;
    }
public:
    size_t get_memory_impact()
    {
        size_t mem_offsets = (num_vectors+1) * sizeof(I);
        size_t mem_vals = size * sizeof(T);
        mem_offsets = utils::round_up(mem_offsets, ator->get_align());
        mem_vals = utils::round_up(mem_vals, ator->get_align());
        size_t mem = mem_offsets + mem_vals;
        return mem;
    }
};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_MULTI_VECTOR_DENSE_DATA_NEW_H_ */
