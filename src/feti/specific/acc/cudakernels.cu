// #ifdef HAVE_CUDA
#include "helper_cuda.h"
#include "cudakernels.h"

namespace espreso {

#define MAX_THREADS 512

/* x(p) = b, for dense vectors x and b */
__global__
void IpvecReorderKernel(const int* __restrict__ p, const double* __restrict__ b, double* __restrict__ x, int n) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid < n) {
        // TODO: This can be optimized by replacing x[tid] = b[tid] with cudamemcpy cudaDeviceToDevice, but it happens only if no reordering applied
        x[p ? p [tid] : tid] = b[tid];
    }
}


__global__
void IpvecReorderMrhsKernel(const int* __restrict__ p, const double* __restrict__ b, double* __restrict__ x, int n, int n_rhs) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int offset = (tid / n) * n;

    if(tid < n * n_rhs) {
        // TODO: This can be optimized by replacing x[tid] = b[tid] with cudamemcpy cudaDeviceToDevice, but it happens only if no reordering applied
        x[p ? p [tid % n] + offset : tid] = b[tid];
    }
}


void cuda::IpvecReorder(const int *p, const double* b, double *x, int n, cudaStream_t stream) {
    int blocks = (n + MAX_THREADS - 1) / MAX_THREADS;
    IpvecReorderKernel<<<blocks, MAX_THREADS, 0, stream>>>(p, b, x, n);
    checkCudaErrors(cudaPeekAtLastError());
}

void cuda::IpvecReorderMrhs(const int *p, const double* b, double *x, int n, int n_rhs, cudaStream_t stream) {
    int blocks = (n * n_rhs + MAX_THREADS - 1) / MAX_THREADS;
    IpvecReorderMrhsKernel<<<blocks, MAX_THREADS, 0, stream>>>(p, b, x, n, n_rhs);
    checkCudaErrors(cudaPeekAtLastError());
}


/* x = b(p), for dense vectors x and b */
__global__
void PvecReorderKernel(const int* __restrict__ p, const double* __restrict__ b, double* __restrict__ x, int n) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid < n) {
        // TODO: This can be optimized by replacing x[tid] = b[tid] with cudamemcpy cudaDeviceToDevice, but it happens only if no reordering applied
        x[tid] = b[p ? p [tid] : tid];
    }
}


__global__
void PvecReorderMrhsKernel(const int* __restrict__ p, const double* __restrict__ b, double* __restrict__ x, int n, int n_rhs) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int offset = (tid / n) * n;

    if(tid < n * n_rhs) {
        // TODO: This can be optimized by replacing x[tid] = b[tid] with cudamemcpy cudaDeviceToDevice, but it happens only if no reordering applied
        x[tid] = b[p ? p [tid % n] + offset : tid];
    }
}


void cuda::PvecReorder(const int *p, const double* b, double *x, int n, cudaStream_t stream) {
    int blocks = (n + MAX_THREADS - 1) / MAX_THREADS;
    PvecReorderKernel<<<blocks, MAX_THREADS, 0, stream>>>(p, b, x, n);
    checkCudaErrors(cudaPeekAtLastError());
}


void cuda::PvecReorderMrhs(const int *p, const double* b, double *x, int n, int n_rhs, cudaStream_t stream) {
    int blocks = (n * n_rhs + MAX_THREADS - 1) / MAX_THREADS;
    PvecReorderMrhsKernel<<<blocks, MAX_THREADS, 0, stream>>>(p, b, x, n, n_rhs);
    checkCudaErrors(cudaPeekAtLastError());
}

}
// #endif