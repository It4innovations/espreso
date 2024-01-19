
#include <cuda_runtime.h>

#ifdef HAVE_CUDA

namespace espreso {

template<typename T, typename I>
static __global__ void _do_DCmap_scatter(T ** domain_vectors, const I * n_dofs_interfaces, const T * cluster_vector, I const * const * D2Cs)
{
    // one block per domain

    I d = blockIdx.x;
    I n_dofs_interface = n_dofs_interfaces[d];
    T * domain_vector = domain_vectors[d];
    const I * D2C = D2Cs[d];

    for(I dof = threadIdx.x; dof < n_dofs_interface; dof += blockDim.x)
    {
        domain_vector[dof] = cluster_vector[D2C[dof]];
    }
}

template<typename T, typename I>
static __global__ void _do_DCmap_gather(T const * const * domain_vectors, const I * n_dofs_interfaces, T * cluster_vector, I const * const * D2Cs)
{
    // one block per domain

    I d = blockIdx.x;
    I n_dofs_interface = n_dofs_interfaces[d];
    const T * domain_vector = domain_vectors[d];
    const I * D2C = D2Cs[d];

    for(I dof = threadIdx.x; dof < n_dofs_interface; dof += blockDim.x)
    {
        atomicAdd(&cluster_vector[D2C[dof]], domain_vector[dof]);
    }
}

template<typename T, typename I>
void _launch_do_DCmap_scatter(dim3 grid_dim, dim3 block_dim, size_t shmem_size, cudaStream_t stream, T ** domain_vectors, const I * n_dofs_interfaces, const T * cluster_vector, I const * const * D2Cs)
{
    _do_DCmap_scatter<T,I><<< grid_dim, block_dim, shmem_size, stream >>>(domain_vectors, n_dofs_interfaces, cluster_vector, D2Cs);
}

template<typename T, typename I>
void _launch_do_DCmap_gather(dim3 grid_dim, dim3 block_dim, size_t shmem_size, cudaStream_t stream, T const * const * domain_vectors, const I * n_dofs_interfaces, T * cluster_vector, I const * const * D2Cs)
{
    _do_DCmap_gather<T,I><<< grid_dim, block_dim, shmem_size, stream >>>(domain_vectors, n_dofs_interfaces, cluster_vector, D2Cs);
}

#define INSTANTIATE(T,I) template void _launch_do_DCmap_scatter<T,I>(dim3 grid_dim, dim3 block_dim, size_t shmem_size, cudaStream_t stream, T ** domain_vectors, const I * n_dofs_interfaces, const T * cluster_vector, I const * const * D2Cs)
INSTANTIATE(float,  int32_t);
//INSTANTIATE(double, int32_t);
INSTANTIATE(float,  int64_t);
//INSTANTIATE(double, int64_t);
#undef INSTANTIATE

#define INSTANTIATE(T,I) template void _launch_do_DCmap_gather<T,I>(dim3 grid_dim, dim3 block_dim, size_t shmem_size, cudaStream_t stream, T const * const * domain_vectors, const I * n_dofs_interfaces, T * cluster_vector, I const * const * D2Cs)
INSTANTIATE(float,  int32_t);
//INSTANTIATE(double, int32_t);
INSTANTIATE(float,  int64_t);
//INSTANTIATE(double, int64_t);
#undef INSTANTIATE

}

#endif // HAVE_CUDA 