
#ifdef HAVE_CUDA

#include "wrappers/cuda/common_cusparse.h"



namespace espreso {
namespace gpu {



handle_cusparse_new::handle_cusparse_new()
{
    internal = std::make_unique<handle_cusparse_new_internal>();
}



handle_cusparse_new::~handle_cusparse_new()
{
    internal.reset();
}



}
}

#endif
