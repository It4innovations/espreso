
#ifndef SRC_GPU_MANAGEMENT_H_
#define SRC_GPU_MANAGEMENT_H_

#include <vector>
#include <functional>
#include <cstring>

#include "math/primitives/vector_dense.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"

namespace espreso {
namespace gpu {
namespace mgm {

    struct _device;
    using device = std::shared_ptr<_device>;

    struct _queue;
    using queue = std::shared_ptr<_queue>;

    device get_device_by_mpi(int mpi_rank, int mpi_size);

    void init_gpu(device & d);

    void set_device(device & d); // global variable representing the device thas is being used

    void queue_create(queue & q);

    void queue_destroy(queue & q);

    void queue_async_barrier(const std::vector<queue> & waitfor, const std::vector<queue> & waitin);

    void queue_wait(queue & q);

    void device_wait();

    size_t get_device_memory_capacity();

    void * memalloc_device(size_t num_bytes);

    void memfree_device(void * ptr);
    
    void memalloc_device_max(void * & memory, size_t & memory_size_B, size_t max_needed);

    void * memalloc_hostpinned(size_t num_bytes);

    void memfree_hostpinned(void * ptr);

    void submit_host_function(queue & q, const std::function<void(void)> & f);

    template<typename T, typename I>
    void copy_submit_h2d(queue & q, T * dst, T const * src, I num_elements);

    template<typename T, typename I>
    void copy_submit_d2h(queue & q, T * dst, T const * src, I num_elements);

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input);

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Vector_Dense<T,I,Ao> & output, const Vector_Dense<T,I,Ai> & input);

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input);

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input);

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_h2d(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern = true, bool copy_vals = true);
    
    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit_d2h(queue & q, Matrix_CSR<T,I,Ao> & output, const Matrix_CSR<T,I,Ai> & input, bool copy_pattern = true, bool copy_vals = true);

    void memset_submit(queue & q, void * ptr, size_t num_bytes, char val);

    inline char change_operation_array_transpose(char op)
    {
        switch(op)
        {
            case 'N': return 'T';
            case 'T': return 'N';
            case 'C': return 'H';
            case 'H': return 'C';
            default: return '_';
        }
    }

    inline char change_operation_conj_transpose(char op)
    {
        switch(op)
        {
            case 'N': return 'H';
            case 'T': return 'C';
            case 'C': return 'T';
            case 'H': return 'N';
            default: return '_';
        }
    }

    inline char change_order(char order)
    {
        if(order == 'R') return 'C';
        if(order == 'C') return 'R';
        return '_';
    }

    class Ad // Allocator device
    {
    public:
        static constexpr bool is_data_host_accessible = false;
        static constexpr bool is_data_device_accessible = true;
        static constexpr bool always_equal = true; // data are allocated from the same globally-set device
    public:
        void * allocate(size_t num_bytes)
        {
            return memalloc_device(num_bytes);
        }
        void deallocate(void * ptr)
        {
            memfree_device(ptr);
        }
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
        static constexpr bool always_equal = true;
    public:
        void * allocate(size_t num_bytes)
        {
            return memalloc_hostpinned(num_bytes);
        }
        void deallocate(void * ptr)
        {
            memfree_hostpinned(ptr);
        }
        template<typename T>
        T * allocate(size_t count)
        {
            return reinterpret_cast<T*>(allocate(count * sizeof(T)));
        }
    };

    template<typename T, typename I, typename A>
    void print_matrix_dense(queue & q, Matrix_Dense<T,I,A> & matrix, const char * name = "")
    {
        if constexpr(A::is_data_host_accessible)
        {
            printf("Dense matrix %s, size %lldx%lld, ld %lld\n", name, (long long)matrix.nrows, (long long)matrix.ncols, (long long)matrix.get_ld());
            for(I r = 0; r < matrix.nrows; r++)
            {
                for(I c = 0; c < matrix.ncols; c++)
                {
                    if constexpr(std::is_floating_point_v<T>)
                    {
                        double v = (double)matrix.vals[r * matrix.get_ld() + c];
                        char str[100];
                        snprintf(str, sizeof(str), "%+11.3e", v);
                        if(strstr(str, "nan") != nullptr) printf("   nan      ");
                        else if(strstr(str, "inf") != nullptr) printf("  %cinf      ", v > 0 ? '+' : '-');
                        else if(v == 0) printf("   0        ");
                        else printf(" %+11.3e", v);
                    }
                    if constexpr(std::is_integral_v<T>) printf(" %+11lld", (long long)matrix.vals[r * matrix.get_ld() + c]);
                }
                printf("\n");
            }
            fflush(stdout);
        }
        else if constexpr(A::is_data_device_accessible)
        {
            Matrix_Dense<T,I,Ah> matrix_host;
            matrix_host.resize(matrix);
            copy_submit_d2h(q, matrix_host, matrix);
            queue_wait(q);
            print_matrix_dense(q, matrix_host, name);
        }
        else
        {
            static_assert(true, "weird matrix with inaccessible data");
        }
    }

    template<typename T, typename I, typename A>
    void print_matrix_csr_as_dense(queue & q, Matrix_CSR<T,I,A> & matrix, const char * name = "")
    {
        if constexpr(A::is_data_host_accessible)
        {
            printf("CSR matrix %s, size %lldx%lld, nnz %lld\n", name, (long long)matrix.nrows, (long long)matrix.ncols, (long long)matrix.nnz);
            for(I r = 0; r < matrix.nrows; r++)
            {
                I start = matrix.rows[r];
                I end = matrix.rows[r+1];
                I curr_col = 0;
                for(I i = start; i <= end; i++)
                {
                    I col = matrix.ncols;
                    if(i < end) col = matrix.cols[i];
                    for(I c = curr_col; c < col; c++) printf("    .       ");
                    if(i < end)
                    {
                        double val = matrix.vals[i];
                        if(val == 0) printf("    0       ");
                        else printf(" %+11.3e", val);
                        curr_col = col + 1;
                    }
                }
                printf("\n");
            }
            fflush(stdout);
        }
        else if constexpr(A::is_data_device_accessible)
        {
            Matrix_CSR<T,I,Ah> matrix_host;
            matrix_host.resize(matrix);
            copy_submit_d2h(q, matrix_host, matrix);
            queue_wait(q);
            print_matrix_csr_as_dense(q, matrix_host, name);
        }
        else
        {
            static_assert(true, "weird matrix with inaccessible data");
        }
    }

}
}
}

#endif /* SRC_GPU_MANAGEMENT_H_ */
