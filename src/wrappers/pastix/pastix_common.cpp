
#ifdef HAVE_PASTIX

#include "wrappers/pastix/pastix_common.h"

static int total_pastix_gpu_instances = 0;

void espreso::math::operations::check_pastix_instances(bool use_gpu, bool created)
{
    if(use_gpu && created) {
        #pragma omp critical(pastix_gpu_instances)
        {
            total_pastix_gpu_instances++;
            if(total_pastix_gpu_instances > 1) {
                // limitation with global structures inside starpu and parsec and probably buggy init/destroy behavior
                eslog::error("only one gpu pastix instance can be created in a program\n");
            }
        }
    }
}

#endif
