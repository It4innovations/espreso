
#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"
#include <sstream>
#include <numeric>
#include <algorithm>

namespace espreso {
namespace gpu {
namespace mgm {
    int _internal_get_local_gpu_idx(int device_count)
    {
        MPI_Comm comm_shared;
        MPI_Comm_split_type(info::mpi::comm, MPI_COMM_TYPE_SHARED, info::mpi::rank, MPI_INFO_NULL, &comm_shared);
        int local_mpi_rank, local_mpi_size;
        MPI_Comm_rank(comm_shared, &local_mpi_rank);
        MPI_Comm_size(comm_shared, &local_mpi_size);
        MPI_Comm_free(&comm_shared);

        if(device_count != local_mpi_size) {
            eslog::error("inconsistent mpi-gpu setting. ranks_per_node=%d != gpus_per_node=%d\n", local_mpi_size, device_count);
        }

        const char * env_map_localrank_gpu = utils::getEnv("ESPRESO_MAP_LOCALRANK_TO_GPU");
        std::vector<int> map;
        if(env_map_localrank_gpu) {
            std::istringstream is(env_map_localrank_gpu);
            while(true)
            {
                int gpu_idx;
                is >> gpu_idx;
                if(is.fail()) {
                    break;
                }
                else {
                    map.push_back(gpu_idx);
                }
            }

            if((int)map.size() != local_mpi_size) {
                eslog::error("incorrect ESPRESO_MAP_LOCALRANK_TO_GPU, contains %zu indices, but there are %d gpus on node\n", map.size(), device_count);
            }
            if(std::any_of(map.begin(), map.end(), [device_count](int idx){ return idx < 0 || idx >= device_count; })) {
                eslog::error("incorrect ESPRESO_MAP_LOCALRANK_TO_GPU, out-of-range values\n");
            }
            int count = local_mpi_size;
            // test if it contains all integers 0..(N-1)
            if(std::accumulate(map.begin(), map.end(), 0) != count * (count-1) / 2) {
                eslog::error("incorrect ESPRESO_MAP_LOCALRANK_TO_GPU, failed sum-of-elements test\n");
            }
            if(std::transform_reduce(map.begin(), map.end(), 0, std::plus<int>(), [](int a){return a*a;}) != (count-1) * count * (2*count-1) / 6) {
                eslog::error("incorrect ESPRESO_MAP_LOCALRANK_TO_GPU, failed sum-of-squares test\n");
            }
        }
        else {
            for(int i = 0; i < local_mpi_size; i++) {
                map.push_back(i);
            }
        }

        return map[local_mpi_rank];
    }
}
}
}



#ifndef ESPRESO_USE_WRAPPER_GPU_CUDA
#ifndef ESPRESO_USE_WRAPPER_GPU_ROCM
#ifndef ESPRESO_USE_WRAPPER_GPU_ONEAPI

#include "gpu_management.h"
#include "basis/utilities/cbmb_allocator.h"

namespace espreso {
namespace gpu {
namespace mgm {

    gpu_wrapper_impl get_implementation()
    {
        return gpu_wrapper_impl::NONE;
    }

    bool is_available() { return false; }

    struct _device {};

    struct _queue {};

    device get_device_by_mpi(int /*mpi_rank*/, int /*mpi_size*/) { eslog::error("calling empty GPU wrapper.\n"); }

    void init_gpu(device & /*d*/) { eslog::error("calling empty GPU wrapper.\n"); }

    void set_device(device & /*d*/) { eslog::error("calling empty GPU wrapper.\n"); }

    void queue_create(queue & /*q*/) { eslog::error("calling empty GPU wrapper.\n"); }

    void queue_destroy(queue & /*q*/) {}

    void queue_async_barrier(const std::vector<queue> & /*waitfor*/, const std::vector<queue> & /*waitin*/) {}

    void queue_wait(queue & /*q*/) {}

    void device_wait() {}

    size_t get_device_memory_capacity() { return 0; }

    size_t get_device_memory_free() { return 0; }

    size_t get_natural_pitch_align() { return 1; }

    void * memalloc_device(size_t /*num_bytes*/) { return nullptr; }

    void * memalloc_device_2d(size_t /*num_chunks*/, size_t /*bytes_per_chunk*/, size_t & /*pitch*/) { return nullptr; }

    void memfree_device(void * /*ptr*/) {}

    void memalloc_device_max(void * & /*memory*/, size_t & /*memory_size_B*/, size_t /*max_needed*/) {}

    void * memalloc_hostpinned(size_t /*num_bytes*/) { return nullptr; }

    void memfree_hostpinned(void * /*ptr*/) {}

    void submit_host_function(queue & /*q*/, const std::function<void(void)> & /*f*/) {}

    template<typename T>
    void copy_submit(queue & /*q*/, T * /*dst*/, T const * /*src*/, size_t /*num_elements*/) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit(queue & /*q*/, Vector_Dense<T,I,Ao> & /*output*/, const Vector_Dense<T,I,Ai> & /*input*/) {}

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit(queue & /*q*/, Matrix_Dense<T,I,Ao> & /*output*/, const Matrix_Dense<T,I,Ai> & /*input*/) {}
    
    template<typename T, typename I, typename Ao, typename Ai>
    void copy_submit(queue & /*q*/, Matrix_CSR<T,I,Ao> & /*output*/, const Matrix_CSR<T,I,Ai> & /*input*/, bool /*copy_pattern*/, bool /*copy_vals*/) {}

    template<typename T>
    void copy_submit(queue & /*q*/, PermutationView_new<T> & /*src*/, PermutationView_new<T> & /*dst*/) {}

    template<typename T>
    void copy_submit(queue & /*q*/, VectorDenseView_new<T> & /*src*/, VectorDenseView_new<T> & /*dst*/) {}

    template<typename T, typename I>
    void copy_submit(queue & /*q*/, MultiVectorDenseView_new<T,I> & /*src*/, MultiVectorDenseView_new<T,I> & /*dst*/, bool /*copy_pattern*/, bool /*copy_vals*/) {}

    template<typename T>
    void copy_submit(queue & /*q*/, MatrixDenseView_new<T> & /*src*/, MatrixDenseView_new<T> & /*dst*/) {}

    template<typename T, typename I>
    void copy_submit(queue & /*q*/, MatrixCsxView_new<T,I> & /*src*/, MatrixCsxView_new<T,I> & /*dst*/, bool /*copy_pattern*/, bool /*copy_vals*/) {}

    void memset_submit(queue & /*q*/, void * /*ptr*/, size_t /*num_bytes*/, char /*val*/) {}

}
}
}

#include "gpu/gpu_management.inst.hpp"

#endif
#endif
#endif
