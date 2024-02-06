
#ifdef HAVE_CUDA

#include "w.cuda.gpu_management.hpp"

namespace espreso {
namespace gpu {
namespace mgm {

device get_device_by_mpi(int mpi_rank, int mpi_size)
{
    device d;
    d.inner = new _device();
    *d.inner = _get_device_by_mpi(mpi_rank, mpi_size);
    return d;
}

void init_gpu(device & d)
{
    d.inner = new _device();
    init_gpu(*d.inner);
}

void set_device(const device & d)
{
    set_device(*d.inner);
}

void queue_create(queue & q, device & d)
{
    q.inner = new _queue();
    queue_create(*q.inner, *d.inner);
}

void queue_destroy(queue & q)
{
    queue_destroy(*q.inner);
    delete q.inner;
}

//void queue_async_barrier(const std::vector<queue> & waitfor, const std::vector<queue> & waitin)
//{
//    queue_async_barrier(waitfor, waitin);
//}

void queue_wait(queue & q)
{
    queue_wait(*q.inner);
}

void device_wait(device & d)
{
    device_wait(*d.inner);
}

size_t get_device_memory_capacity(device & d)
{
    return get_device_memory_capacity(*d.inner);
}

void * memalloc_device(device & d, size_t num_bytes)
{
    return memalloc_device(*d.inner, num_bytes);
}

void memfree_device(device & d, void * ptr)
{
    memfree_device(*d.inner, ptr);
}

void memalloc_device_max(device & d, void * & memory, size_t & memory_size_B, size_t max_needed)
{
    memalloc_device_max(*d.inner, memory, memory_size_B, max_needed);
}

void submit_host_function(queue & q, const std::function<void(void)> & c)
{
    submit_host_function(*q.inner, c);
}

template<> void copy_submit_d2h<double , int                              >(queue & q, double * dst, const double * src, int num_elements) { copy_submit_d2h(*q.inner, dst, src, num_elements); }
template<> void copy_submit_d2h<double , int, mgm::Ah      , mgm::Ad      >(queue & q, Vector_Dense<double , int, mgm::Ah      > & output, const Vector_Dense<double , int, mgm::Ad      > & input) { copy_submit_d2h(*q.inner, output, input); }
template<> void copy_submit_d2h<double , int, cpu_allocator, mgm::Ad      >(queue & q, Vector_Dense<double , int, cpu_allocator> & output, const Vector_Dense<double , int, mgm::Ad      > & input) { copy_submit_d2h(*q.inner, output, input); }

template<> void copy_submit_h2d<int    , int                              >(queue & q, int    * dst, const int    * src, int num_elements) { copy_submit_h2d(*q.inner, dst, src, num_elements); }
template<> void copy_submit_h2d<double , int                              >(queue & q, double * dst, const double * src, int num_elements) { copy_submit_h2d(*q.inner, dst, src, num_elements); }
template<> void copy_submit_h2d<double , int, mgm::Ad      , mgm::Ah      >(queue & q, Vector_Dense<double , int, mgm::Ad      > & output, const Vector_Dense<double , int, mgm::Ah      > & input) { copy_submit_h2d(*q.inner, output, input); }
template<> void copy_submit_h2d<int    , int, mgm::Ad      , mgm::Ah      >(queue & q, Vector_Dense<int    , int, mgm::Ad      > & output, const Vector_Dense<int    , int, mgm::Ah      > & input) { copy_submit_h2d(*q.inner, output, input); }
template<> void copy_submit_h2d<double , int, mgm::Ad      , cpu_allocator>(queue & q, Vector_Dense<double , int, mgm::Ad      > & output, const Vector_Dense<double , int, cpu_allocator> & input) { copy_submit_h2d(*q.inner, output, input); }
template<> void copy_submit_h2d<int*   , int, mgm::Ad      , mgm::Ah      >(queue & q, Vector_Dense<int*   , int, mgm::Ad      > & output, const Vector_Dense<int*   , int, mgm::Ah      > & input) { copy_submit_h2d(*q.inner, output, input); }
template<> void copy_submit_h2d<double*, int, mgm::Ad      , mgm::Ah      >(queue & q, Vector_Dense<double*, int, mgm::Ad      > & output, const Vector_Dense<double*, int, mgm::Ah      > & input) { copy_submit_h2d(*q.inner, output, input); }
template<> void copy_submit_h2d<double , int, mgm::Ad      , mgm::Ah      >(queue & q, Matrix_Dense<double , int, mgm::Ad      > & output, const Matrix_Dense<double , int, mgm::Ah      > & input) { copy_submit_h2d(*q.inner, output, input); }
template<> void copy_submit_h2d<double , int, mgm::Ad      , mgm::Ah      >(queue & q, Matrix_CSR  <double , int, mgm::Ad      > & output, const Matrix_CSR  <double , int, mgm::Ah      > & input, bool copy_pattern, bool copy_vals) { copy_submit_h2d(*q.inner, output, input, copy_pattern, copy_vals); }

template<>
void memset_submit<>(queue & q, void * ptr, size_t num_bytes, char val)
{
    memset_submit<double, size_t>(*q.inner, ptr, num_bytes, val);
}

}
}
}

#endif
