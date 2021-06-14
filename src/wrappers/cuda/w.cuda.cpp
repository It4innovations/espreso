
#include "w.cuda.h"

#include "esinfo/eslog.hpp"
#include "wrappers/mpi/communication.h"

int espreso::cuda::devices = 0;
size_t espreso::cuda::availMemory = 0;
size_t espreso::cuda::totalMemory = 0;

#ifdef HAVE_CUDA
#include "cuda_runtime_api.h"
#include "helper_cuda.h"

namespace espreso {

void cuda::fillDeviceInfo()
{
	cudaError_t cuda_status = cudaGetDeviceCount(&cuda::devices);
	if(cuda_status != cudaSuccess) {
		espreso::cuda::devices = 0;
		eslog::error("ERROR: No CUDA device available (--solver=cuda or --with-cuda used).\n Underlying CUDA error at %s:%d code=%d(%s) \"%s\"\n Descr: %s\n", __FILE__, __LINE__,
            static_cast<unsigned int>(cuda_status), _cudaGetErrorEnum(cuda_status), "cudaGetDeviceCount(&cuda::devices)", cudaGetErrorString((cudaError_t)cuda_status));
	}

	if (MPITools::node->rank == 0) {
		for (int i = 0; i < cuda::devices; i++) {
			cudaDeviceProp prop;
			checkCudaErrors(cudaGetDeviceProperties(&prop, i));
			checkCudaErrors(cudaMemGetInfo(&cuda::availMemory, &cuda::totalMemory));
			cuda::availMemory /= 1024 * 1024;
			cuda::totalMemory /= 1024 * 1024;
		}
	}
	size_t tmp[3] = { (size_t)cuda::devices, cuda::availMemory, cuda::totalMemory };
	Communication::allReduce(tmp, NULL, 3, MPITools::getType<size_t>().mpitype, MPI_SUM);

	eslog::info(" == NUMBER OF ACCELERATORS %*d == \n", 64, cuda::devices);
	eslog::info(" == ACCELERATORS MEMORY :: AVAILABLE %*d [MB] == \n", 49, cuda::availMemory);
	eslog::info(" == ACCELERATORS MEMORY :: TOTAL %*d [MB] == \n", 53, cuda::totalMemory);
	eslog::info(" ============================================================================================= \n");
}

void cuda::SetDevice(esint device_id) {
	checkCudaErrors(cudaSetDevice(device_id));
}

void cuda::Malloc(void **dev_ptr, size_t size) {
	checkCudaErrors(cudaMalloc(dev_ptr, size));
}

}
#endif
