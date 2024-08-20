
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

    enum struct gpu_wrapper_impl {
        NONE,
        CUDA,
        ROCM,
        ONEAPI
    };

    gpu_wrapper_impl get_implementation();

    inline bool is_linked()
    {
        return (get_implementation() != gpu_wrapper_impl::NONE);
    }

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

    void event_record(event & e, queue & q);

    void event_wait(event & e);

    size_t get_device_memory_capacity();

    void * memalloc_device(size_t num_bytes);

    void memfree_device(void * ptr);
    
    void memalloc_device_max(void * & memory, size_t & memory_size_B, size_t max_needed);

    void * memalloc_hostpinned(size_t num_bytes);

    void memfree_hostpinned(void * ptr);

    void submit_host_function(queue & q, const std::function<void(void)> & f);

    template<typename T>
    void copy_submit_h2d(queue & q, T * dst, T const * src, size_t num_elements);

    template<typename T>
    void copy_submit_d2h(queue & q, T * dst, T const * src, size_t num_elements);

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

    inline char operation_combine(char op1, char op2)
    {
        bool trans1 = (op1 == 'T' || op1 == 'H');
        bool trans2 = (op2 == 'T' || op2 == 'H');
        bool conj1 = (op1 == 'C' || op1 == 'H');
        bool conj2 = (op2 == 'C' || op2 == 'H');
        bool trans = (trans1 != trans2);
        bool conj = (conj1 != conj2);
        if(!trans && !conj) return 'N';
        if(!trans &&  conj) return 'C';
        if( trans && !conj) return 'T';
        if( trans &&  conj) return 'H';
        return '_';
    }

    inline char operation_remove_conj(char op)
    {
        switch(op) {
            case 'N':
            case 'C':
                return 'N';
            case 'T':
            case 'H':
                return 'T';
            default:
                return '_';
        }
    }

    inline char order_change(char order)
    {
        if(order == 'R') return 'C';
        if(order == 'C') return 'R';
        return '_';
    }

    inline char side_change(char side)
    {
        if(side == 'L') return 'R';
        if(side == 'R') return 'L';
        return '_';
    }

    inline char fill_change(char fill)
    {
        if(fill == 'U') return 'L';
        if(fill == 'L') return 'U';
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
    void print_vector(queue & q, Vector_Dense<T,I,A> & vec, const char * name = "")
    {
        if constexpr(A::is_data_host_accessible) {
            printf("Vector %s, size %lld\n", name, (long long)vec.size);
            for(I i = 0; i < vec.size; i++) {
                if constexpr(std::is_floating_point_v<T>) {
                    double v = (double)vec.vals[i];
                    char str[100];
                    snprintf(str, sizeof(str), "%+11.3e", v);
                    if(strstr(str, "nan") != nullptr) printf("   nan      ");
                    else if(strstr(str, "inf") != nullptr) printf("  %cinf      ", v > 0 ? '+' : '-');
                    else if(v == 0) printf("   0        ");
                    else printf(" %+11.3e", v);
                }
                if constexpr(std::is_integral_v<T>) printf(" %+11lld", (long long)vec.vals[i]);
            }
            printf("\n");
            fflush(stdout);
        }
        else if constexpr(A::is_data_device_accessible) {
            Vector_Dense<T,I,Ah> vec_host;
            vec_host.resize(vec);
            copy_submit_d2h(q, vec_host, vec);
            queue_wait(q);
            print_vector(q, vec_host, name);
        }
        else {
            static_assert(true, "weird vector with inaccessible data");
        }
    }

    template<typename T, typename I, typename A>
    void print_matrix_dense(queue & q, Matrix_Dense<T,I,A> & matrix, const char * name = "")
    {
        if constexpr(A::is_data_host_accessible) {
            printf("Dense matrix %s, size %lldx%lld, ld %lld\n", name, (long long)matrix.nrows, (long long)matrix.ncols, (long long)matrix.get_ld());
            for(I r = 0; r < matrix.nrows; r++) {
                for(I c = 0; c < matrix.ncols; c++) {
                    if constexpr(std::is_floating_point_v<T>) {
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
        else if constexpr(A::is_data_device_accessible) {
            Matrix_Dense<T,I,Ah> matrix_host;
            matrix_host.resize(matrix);
            copy_submit_d2h(q, matrix_host, matrix);
            queue_wait(q);
            print_matrix_dense(q, matrix_host, name);
        }
        else {
            static_assert(true, "weird matrix with inaccessible data");
        }
    }

    template<typename T, typename I, typename A>
    void print_matrix_csr_as_dense(queue & q, Matrix_CSR<T,I,A> & matrix, const char * name = "")
    {
        if constexpr(A::is_data_host_accessible) {
            printf("CSR matrix %s, size %lldx%lld, nnz %lld\n", name, (long long)matrix.nrows, (long long)matrix.ncols, (long long)matrix.nnz);
            for(I r = 0; r < matrix.nrows; r++) {
                I start = matrix.rows[r];
                I end = matrix.rows[r+1];
                I curr_col = 0;
                for(I i = start; i <= end; i++) {
                    I col = matrix.ncols;
                    if(i < end) col = matrix.cols[i];
                    for(I c = curr_col; c < col; c++) {
                        printf("    .       ");
                    }
                    if(i < end) {
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
        else if constexpr(A::is_data_device_accessible) {
            Matrix_CSR<T,I,Ah> matrix_host;
            matrix_host.resize(matrix);
            copy_submit_d2h(q, matrix_host, matrix);
            queue_wait(q);
            print_matrix_csr_as_dense(q, matrix_host, name);
        }
        else {
            static_assert(true, "weird matrix with inaccessible data");
        }
    }

    template<typename T, typename I, typename A>
    void print_matrix_csr_arrays(queue & q, Matrix_CSR<T,I,A> & matrix, const char * name = "")
    {
        if constexpr(A::is_data_host_accessible) {
            printf("CSR matrix %s, size %lldx%lld, nnz %lld\n", name, (long long)matrix.nrows, (long long)matrix.ncols, (long long)matrix.nnz);
            printf("row ptrs: ");
            for(I r = 0; r <= matrix.nrows; r++) printf("%lld ", (long long)matrix.rows[r]);
            printf("\n");
            printf("col idxs: ");
            for(I i = 0; i < matrix.nnz; i++) printf("%lld ", (long long)matrix.cols[i]);
            printf("\n");
            printf("vals:     ");
            for(I i = 0; i < matrix.nnz; i++) printf("%+.3e ", (double)matrix.vals[i]);
            printf("\n");
            fflush(stdout);
        }
        else if constexpr(A::is_data_device_accessible) {
            Matrix_CSR<T,I,Ah> matrix_host;
            matrix_host.resize(matrix);
            copy_submit_d2h(q, matrix_host, matrix);
            queue_wait(q);
            print_matrix_csr_arrays(q, matrix_host, name);
        }
        else {
            static_assert(true, "weird matrix with inaccessible data");
        }
    }

}
}
}

#endif /* SRC_GPU_MANAGEMENT_H_ */
