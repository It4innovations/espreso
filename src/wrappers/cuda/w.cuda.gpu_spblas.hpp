
#ifndef SRC_WRAPPERS_CUDA_W_CUDA_GPU_SPBLAS_HPP_
#define SRC_WRAPPERS_CUDA_W_CUDA_GPU_SPBLAS_HPP_

#ifdef HAVE_CUDA

#ifdef USE_CUSPARSE_LEGACY
#include "w.cuda.gpu_spblas_legacy.hpp"
#else
#include "w.cuda.gpu_spblas_modern.hpp"
#endif

#endif

#endif /* SRC_WRAPPERS_CUDA_W_CUDA_GPU_SPBLAS_HPP_ */
