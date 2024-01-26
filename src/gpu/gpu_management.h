
#ifndef SRC_GPU_MANAGEMENT_H_
#define SRC_GPU_MANAGEMENT_H_

#include <vector>

#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"

namespace espreso {
namespace gpu {
namespace mgm {

    struct device;

    struct queue;


    class Ad // Allocator device
    {
    public:
        static constexpr bool is_data_host_accessible = false;
        static constexpr bool is_data_device_accessible = true;
    public:
        void * allocate(size_t num_bytes);
        template<typename T>
        T * allocate(size_t count);
        template<typename T>
        void deallocate(T * ptr);
    };

    class Ah // Allocator host
    {
    public:
        static constexpr bool is_data_host_accessible = true;
        static constexpr bool is_data_device_accessible = false;
    public:
        void * allocate(size_t num_bytes);
        template<typename T>
        T * allocate(size_t count);
        template<typename T>
        void deallocate(T * ptr);
    };

    static device get_device_by_mpi(int mpi_rank);

    static void set_device(const device & d);

    static void queue_create(queue & q, device & d);

    static void queue_destroy(queue & q);

    static void queue_async_barrier(const std::vector<queue> & waitfor, const std::vector<queue> & waitin);

    static void queue_wait(queue & q);

    static void device_wait(device & d);

    static size_t get_device_memory_capacity(device & d);

    static void * memalloc_device(device & d, size_t num_bytes);

    template<typename T>
    static void memfree_device(device & d, T * ptr);
    
    static void memalloc_device_max(device & d, void * & memory, size_t & memory_size_B, size_t max_needed);

    template<typename C>
    static void submit_host_function(queue & q, const C & c);

    template<typename T, typename I>
    static void copy_submit_h2d(queue & q, T * dst, const T * src, I num_elements);

    template<typename T, typename I>
    static void copy_submit_d2h(queue & q, T * dst, const T * src, I num_elements);

    template<typename T, typename I, typename Ao, typename Ai>
    static void copy_matrix_submit_d2h(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input);

    template<typename T, typename I, typename Ao, typename Ai>
    static void copy_matrix_submit_h2d(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input);
    
    template<typename T, typename I, typename Ao, typename Ai>
    static void copy_matrix_submit_d2h(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern = true, bool copy_vals = true);

    template<typename T, typename I, typename Ao, typename Ai>
    static void copy_matrix_submit_h2d(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern = true, bool copy_vals = true);

    template<typename I>
    static void memset_submit(queue & q, void * ptr, I num_bytes, char val);

    static char change_operation(char op)
    {
        if(op == 'N') return 'T';
        if(op == 'T') return 'N';
        return '_';
    }

    static char change_order(char order)
    {
        if(order == 'R') return 'C';
        if(order == 'C') return 'R';
        return '_';
    }

}
}
}

#ifdef HAVE_CUDA
#include "wrappers/cuda/w.cuda.gpu_management.hpp"
#endif

#endif /* SRC_GPU_MANAGEMENT_H_ */
