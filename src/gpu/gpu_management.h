
#ifndef SRC_GPU_MANAGEMENT_H_
#define SRC_GPU_MANAGEMENT_H_

#include <vector>

#include "math/primitives/vector_dense.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"

#include <functional>

namespace espreso {
namespace gpu {
namespace mgm {

    struct _device;
    struct device { _device *inner; };

    struct _queue;
    struct queue { _queue *inner; };

    class Ad // Allocator device
    {
    public:
        static constexpr bool is_data_host_accessible = false;
        static constexpr bool is_data_device_accessible = true;
    public:
        void * allocate(size_t num_bytes);
        void deallocate(void * ptr);

        template<typename T>
        T * allocate(size_t count)
        {
            return reinterpret_cast<T*>(allocate(count * sizeof(T)));
        }
    };

    class Ah // Allocator host
    {
    public:
        static constexpr bool is_data_host_accessible = true;
        static constexpr bool is_data_device_accessible = false;
    public:
        void * allocate(size_t num_bytes);
        void deallocate(void * ptr);

        template<typename T>
        T * allocate(size_t count)
        {
            return reinterpret_cast<T*>(allocate(count * sizeof(T)));
        }
    };

    device get_device_by_mpi(int mpi_rank, int mpi_size);

    void init_gpu(device & d);

    void set_device(const device & d);

    void queue_create(queue & q, device & d);

    void queue_destroy(queue & q);

    void queue_async_barrier(const std::vector<queue> & waitfor, const std::vector<queue> & waitin);

    void queue_wait(queue & q);

    void device_wait(device & d);

    size_t get_device_memory_capacity(device & d);

    void * memalloc_device(device & d, size_t num_bytes);

    void memfree_device(device & d, void * ptr);
    
    void memalloc_device_max(device & d, void * & memory, size_t & memory_size_B, size_t max_needed);

    void * memalloc_hostpinned(device & d, size_t num_bytes);

    void memfree_hostpinned(device & d, void * ptr);

    void submit_host_function(queue & q, const std::function<void(void)> & c);

    template<typename T, typename I>
    void copy_submit_h2d(queue & q, T * dst, const T * src, I num_elements);

    template<typename T, typename I>
    void copy_submit_d2h(queue & q, T * dst, const T * src, I num_elements);

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input);

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input);

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input);

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input);
    
    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern = true, bool copy_vals = true);

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern = true, bool copy_vals = true);

    template<typename I>
    void memset_submit(queue & q, void * ptr, I num_bytes, char val);

    inline char change_operation(char op)
    {
        if(op == 'N') return 'T';
        if(op == 'T') return 'N';
        return '_';
    }

    inline char change_order(char order)
    {
        if(order == 'R') return 'C';
        if(order == 'C') return 'R';
        return '_';
    }
}
}
}

#endif /* SRC_GPU_MANAGEMENT_H_ */
