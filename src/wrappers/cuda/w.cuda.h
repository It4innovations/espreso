
#ifndef SRC_WRAPPERS_CUDA_W_CUDA_H_
#define SRC_WRAPPERS_CUDA_W_CUDA_H_

#include <cstddef>

namespace espreso {
namespace cuda {

extern int devices;
extern size_t availMemory;
extern size_t totalMemory;

void fillDeviceInfo();

/**
 * \brief Set device to be used for GPU executions
 *
 * Sets \p device as the current device for the calling host thread.
 * Valid device id's are 0 to (::cudaGetDeviceCount() - 1).
 *
 * Any device memory subsequently allocated from this host thread
 * using ::cudaMalloc(), ::cudaMallocPitch() or ::cudaMallocArray()
 * will be physically resident on \p device.  Any host memory allocated
 * from this host thread using ::cudaMallocHost() or ::cudaHostAlloc()
 * or ::cudaHostRegister() will have its lifetime associated  with
 * \p device.  Any streams or events created from this host thread will
 * be associated with \p device.  Any kernels launched from this host
 * thread using the <<<>>> operator or ::cudaLaunchKernel() will be executed
 * on \p device.
 *
 * This call may be made from any host thread, to any device, and at
 * any time.  This function will do no synchronization with the previous
 * or new device, and should be considered a very low overhead call.
 *
 * \param device - Device on which the active host thread should execute the
 * device code.
 *
 * \return
 * ::cudaSuccess,
 * ::cudaErrorInvalidDevice,
 * ::cudaErrorDeviceAlreadyInUse
 * \notefnerr
 * \note_init_rt
 * \note_callback
 *
 * \sa ::cudaGetDeviceCount, ::cudaGetDevice, ::cudaGetDeviceProperties,
 * ::cudaChooseDevice,
 * ::cuCtxSetCurrent
 */
void SetDevice(esint device_id);

/**
 * \brief Allocate memory on the device
 *
 * Allocates \p size bytes of linear memory on the device and returns in
 * \p *devPtr a pointer to the allocated memory. The allocated memory is
 * suitably aligned for any kind of variable. The memory is not cleared.
 * ::cudaMalloc() returns ::cudaErrorMemoryAllocation in case of failure.
 *
 * The device version of ::cudaFree cannot be used with a \p *devPtr
 * allocated using the host API, and vice versa.
 *
 * \param devPtr - Pointer to allocated device memory
 * \param size   - Requested allocation size in bytes
 *
 * \return
 * ::cudaSuccess,
 * ::cudaErrorInvalidValue,
 * ::cudaErrorMemoryAllocation
 * \notefnerr
 * \note_init_rt
 * \note_callback
 *
 * \sa ::cudaMallocPitch, ::cudaFree, ::cudaMallocArray, ::cudaFreeArray,
 * ::cudaMalloc3D, ::cudaMalloc3DArray,
 * \ref ::cudaMallocHost(void**, size_t) "cudaMallocHost (C API)",
 * ::cudaFreeHost, ::cudaHostAlloc,
 * ::cuMemAlloc
 */
void Malloc(void **devPtr, size_t size);

}
}

#endif /* SRC_WRAPPERS_CUDA_W_CUDA_H_ */

