
#include "wrappers/pastix/pastix_common.h"

static int total_pastix_gpu_instances = 0;
static int pastix_concurrent_instances = 0;

void espreso::math::operations::check_pastix_instances(bool use_gpu, bool created)
{
    #pragma omp critical(pastix_concurrent_instances)
    {
        if(created) {
            pastix_concurrent_instances++;
            if(pastix_concurrent_instances > 1) {
                // sometimes it works, but sometimes it gives wrong results when run concurrently
                // even critical sections dont solve this
                eslog::error("only one pastix instance can be active at a time\n");
            }
        }
        else {
            pastix_concurrent_instances--;
        }
    }

    if(use_gpu) {
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
