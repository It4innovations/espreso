
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
	checkCudaErrors(cudaGetDeviceCount(&cuda::devices));

	if (MPITools::node->rank == 0) {
		for (int i = 0; i < cuda::devices; i++) {
			cudaDeviceProp prop;
			cudaGetDeviceProperties(&prop, i);
			cudaMemGetInfo(&cuda::availMemory, &cuda::totalMemory);
			cuda::availMemory /= 1024 * 1024;
			cuda::totalMemory /= 1024 * 1024;
		}
	}
	size_t tmp[3] = { (size_t)cuda::devices, cuda::availMemory, cuda::totalMemory };
	Communication::allReduce(tmp, NULL, 3, MPITools::getType(tmp).mpitype, MPI_SUM);

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
