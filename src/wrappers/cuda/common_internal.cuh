
#ifndef SRC_WRAPPERS_CUDA_COMMON_INTERNAL_CUH
#define SRC_WRAPPERS_CUDA_COMMON_INTERNAL_CUH



namespace espreso {
namespace gpu {


template<typename T>
static __device__ void myGenericAddScaled(T * dst, T val, T scalar) { *dst += scalar * val; }
template<typename T>
static __device__ void myGenericAddScaled(std::complex<T> * dst, std::complex<T> val, T scalar)
{
    myGenericAddScaled(&reinterpret_cast<T*>(dst)[0], reinterpret_cast<T*>(&val)[0], scalar);
    myGenericAddScaled(&reinterpret_cast<T*>(dst)[1], reinterpret_cast<T*>(&val)[1], scalar);
}

template<typename T>
static __device__ void myGenericAdd(T * dst, T val) { *dst += val; }
template<typename T>
static __device__ void myGenericAdd(std::complex<T> * dst, std::complex<T> val)
{
    myGenericAdd(&reinterpret_cast<T*>(dst)[0], reinterpret_cast<T*>(&val)[0]);
    myGenericAdd(&reinterpret_cast<T*>(dst)[1], reinterpret_cast<T*>(&val)[1]);
}

template<typename T>
static __device__ void complexAtomicAdd(std::complex<T> * dst, std::complex<T> val)
{
    atomicAdd(&reinterpret_cast<T*>(dst)[0], reinterpret_cast<T*>(&val)[0]);
    atomicAdd(&reinterpret_cast<T*>(dst)[1], reinterpret_cast<T*>(&val)[1]);
}

template<typename T>
static __device__ void myAtomicAdd(T * dst, T val) { atomicAdd(dst, val); }
template<typename T>
static __device__ void myAtomicAdd(std::complex<T> * dst, std::complex<T> val) { complexAtomicAdd(dst, val); }



}
}



#endif /* SRC_WRAPPERS_CUDA_COMMON_INTERNAL_CUH */
