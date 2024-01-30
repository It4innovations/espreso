
#ifdef HAVE_CUDA

#ifdef USE_CUSPARSE_LEGACY
#include "w.cuda.gpu_spblas_legacy.hpp"
#else
#include "w.cuda.gpu_spblas_modern.hpp"
#endif

#endif
