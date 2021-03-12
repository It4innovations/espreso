#ifndef SOLVER_SPECIFIC_ACC_CUDA_KERNELS_H_
#define SOLVER_SPECIFIC_ACC_CUDA_KERNELS_H_

#include <cuda_runtime.h>

namespace espreso {
namespace cuda {

void IpvecReorder(const int *p, const double* b, double *x, int n, cudaStream_t stream);
void IpvecReorderMrhs(const int *p, const double* b, double *x, int n, int n_rhs, cudaStream_t stream);
void PvecReorder(const int *p, const double* b, double *x, int n, cudaStream_t stream);
void PvecReorderMrhs(const int *p, const double* b, double *x, int n, int n_rhs, cudaStream_t stream);

}
}

#endif /* SOLVER_SPECIFIC_ACC_CUDA_KERNELS_H_ */