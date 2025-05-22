
#ifndef SRC_MATH_PRIMITIVES_NEW_VECTOR_DENSE_VIEW_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_VECTOR_DENSE_VIEW_NEW_H_

#include "math/primitives_new/vector_base_new.h"
#include "math/primitives_new/allocator_new_base.h"
#include "math/primitives/vector_dense.h"
#include "basis/utilities/utils.h"

#include <cstring>



namespace espreso {



template<typename T>
struct VectorDenseView_new : public VectorBase_new
{
public: // the user promises not to modify these values (I don't want to implement getters everywhere)
    Allocator_new * ator = nullptr;
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
    void set_view(size_t size_, T * vals_, Allocator_new * ator_)
    {
        if(was_set) eslog::error("can only set yet-uninitialized vector view\n");
        size = size_;
        vals = vals_;
        ator = ator_;
        was_set = true;
    }
public:
    template<typename I, typename A>
    static VectorDenseView_new<T> from_old(const Vector_Dense<T,I,A> & V_old)
    {
        VectorDenseView_new<T> V_new;
        V_new.set_view(V_old.size, V_old.vals, AllocatorDummy_new::get_singleton(A::is_data_host_accessible, A::is_data_device_accessible));
        return V_new;
    }
    template<typename I, typename A>
    static Vector_Dense<T,I,A> to_old(const VectorDenseView_new<T> & V_new)
    {
        if(A::is_data_host_accessible != V_new.ator->is_data_accessible_cpu()) eslog::error("allocator access mismatch on cpu\n");
        if(A::is_data_device_accessible != V_new.ator->is_data_accessible_gpu()) eslog::error("allocator access mismatch on gpu\n");
        Vector_Dense<T,I,A> V_old;
        V_old._allocated.size = V_new.size;
        V_old.size = V_new.size;
        V_old.vals = V_new.vals;
        return V_old;
    }
public:
    void print(const char * name = "")
    {
        if(!ator->is_data_accessible_cpu()) eslog::error("print is supported only for cpu-accessible matrices\n");
        if constexpr(utils::is_real<T>()) {
            eslog::info("Dense vector %s, size %zu\n", name, size);
            eslog::info("vals: ");
            for(size_t i = 0; i < size; i++) {
                if constexpr(std::is_floating_point_v<T>) {
                    double v = (double)vals[i];
                    char str[100];
                    snprintf(str, sizeof(str), "%+11.3e", v);
                    if(strstr(str, "nan") != nullptr) eslog::info("   nan      ");
                    else if(strstr(str, "inf") != nullptr) eslog::info("  %cinf      ", v > 0 ? '+' : '-');
                    else if(v == 0) eslog::info("   0        ");
                    else eslog::info(" %+11.3e", v);
                }
                if constexpr(std::is_integral_v<T>) {
                    long long v = (long long)vals[i];
                    eslog::info(" %+11lld", v);
                }
            }
            eslog::info("\n");
            fflush(stdout);
        }
        if constexpr(utils::is_complex<T>()) {
            eslog::error("vector print not yet supported for complex matrices\n");
        }
    }
};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_VECTOR_DENSE_VIEW_NEW_H_ */
