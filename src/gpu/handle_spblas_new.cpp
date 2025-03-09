
#include "gpu/handle_spblas_new.h"

#include "wrappers/cuda/handle_cusparse_new.h"



namespace espreso {
namespace gpu {



std::unique_ptr<handle_spblas_new> handle_spblas_new::make()
{
    // feel free to make this runtime ifs based on ecf or env
    #ifdef ESPRESO_USE_WRAPPER_GPU_CUDA
        return std::make_unique<handle_cusparse_new>();
    #endif
    eslog::error("wrapper for handle_spblas_new not available");
}



}
}

#endif /* SRC_GPU_HANDLE_SPBLAS_NEW_H */
