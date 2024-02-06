
#include "gpu_management.h"

namespace espreso {
namespace gpu {
namespace mgm {

#ifndef HAVE_CUDA

void* Ad::allocate(size_t num_bytes)
{
    return nullptr;
}

void Ad::deallocate(void *ptr)
{

}

void* Ah::allocate(size_t num_bytes)
{
    return nullptr;
}

void Ah::deallocate(void *ptr)
{

}

device get_device_by_mpi(int mpi_rank, int mpi_size) { return device{}; }

void init_gpu(device & d) {}

void set_device(const device & d) {}

void queue_create(queue & q, device & d) {}

void queue_destroy(queue & q) {}

void queue_async_barrier(const std::vector<queue> & waitfor, const std::vector<queue> & waitin) {}

void queue_wait(queue & q) {}

void device_wait(device & d) {}

size_t get_device_memory_capacity(device & d) { return 0; }

void * memalloc_device(device & d, size_t num_bytes) { return nullptr; }

void memfree_device(device & d, void * ptr) {}

void memalloc_device_max(device & d, void * & memory, size_t & memory_size_B, size_t max_needed) {}

void submit_host_function(queue & q, const std::function<void(void)> & c) {}

template<> void copy_submit_d2h<double , int                              >(queue & q, double * dst, const double * src, int num_elements) {}
template<> void copy_submit_d2h<double , int, mgm::Ah      , mgm::Ad      >(queue & q, Vector_Dense<double , int, mgm::Ah      > & output, const Vector_Dense<double , int, mgm::Ad      > & input) {}
template<> void copy_submit_d2h<double , int, cpu_allocator, mgm::Ad      >(queue & q, Vector_Dense<double , int, cpu_allocator> & output, const Vector_Dense<double , int, mgm::Ad      > & input) {}

template<> void copy_submit_h2d<int    , int                              >(queue & q, int * dst, const int * src, int num_elements) {}
template<> void copy_submit_h2d<double , int                              >(queue & q, double * dst, const double * src, int num_elements) {}
template<> void copy_submit_h2d<double , int, mgm::Ad      , mgm::Ah      >(queue & q, Vector_Dense<double , int, mgm::Ad      > & output, const Vector_Dense<double , int, mgm::Ah      > & input) {}
template<> void copy_submit_h2d<int    , int, mgm::Ad      , mgm::Ah      >(queue & q, Vector_Dense<int    , int, mgm::Ad      > & output, const Vector_Dense<int    , int, mgm::Ah      > & input) {}
template<> void copy_submit_h2d<double , int, mgm::Ad      , cpu_allocator>(queue & q, Vector_Dense<double , int, mgm::Ad      > & output, const Vector_Dense<double , int, cpu_allocator> & input) {}
template<> void copy_submit_h2d<int*   , int, mgm::Ad      , mgm::Ah      >(queue & q, Vector_Dense<int*   , int, mgm::Ad      > & output, const Vector_Dense<int*   , int, mgm::Ah      > & input) {}
template<> void copy_submit_h2d<double*, int, mgm::Ad      , mgm::Ah      >(queue & q, Vector_Dense<double*, int, mgm::Ad      > & output, const Vector_Dense<double*, int, mgm::Ah      > & input) {}
template<> void copy_submit_h2d<double , int, mgm::Ad      , mgm::Ah      >(queue & q, Matrix_Dense<double , int, mgm::Ad      > & output, const Matrix_Dense<double , int, mgm::Ah      > & input) {}
template<> void copy_submit_h2d<double , int, mgm::Ad      , mgm::Ah      >(queue & q, Matrix_CSR  <double , int, mgm::Ad      > & output, const Matrix_CSR  <double , int, mgm::Ah      > & input, bool copy_pattern, bool copy_vals) {}

template<>
void memset_submit<>(queue & q, void * ptr, size_t num_bytes, char val) {}

#endif

}
}
}


