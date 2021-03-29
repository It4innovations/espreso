
#include "w.cuda.h"
#include "esinfo/eslog.h"

#ifndef HAVE_CUDA

namespace espreso {

void cuda::fillDeviceInfo()
{

}

void cuda::SetDevice(esint device_id) {
	eslog::warning("ESPRESO run-time error: cannot call CUDA library (the library is not linked).\n");
}

void cuda::Malloc(void **dev_ptr, size_t size) {
	eslog::warning("ESPRESO run-time error: cannot call CUDA library (the library is not linked).\n");
}

// TODO void haveCUDA() {}

}

#endif
