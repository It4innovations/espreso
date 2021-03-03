#include "w.cuda.h"

#ifdef HAVE_CUDA
#include "cuda_runtime_api.h"
#include "helper_cuda.h"

namespace espreso {

void cuda::SetDevice(esint device_id) {
	checkCudaErrors(cudaSetDevice(device_id));
}

void cuda::Malloc(void **dev_ptr, size_t size) {
	checkCudaErrors(cudaMalloc(dev_ptr, size));
}

}
#endif