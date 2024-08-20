
#ifndef SRC_WRAPPERS_ONEAPI_W_ONEAPI_GPU_MANAGEMENT_H_
#define SRC_WRAPPERS_ONEAPI_W_ONEAPI_GPU_MANAGEMENT_H_

#include <sycl/sycl.hpp>



namespace espreso {
namespace gpu {
namespace mgm {

    struct _device
    {
        sycl::device d;
        sycl::context c;
        std::vector<queue> qs;
        _device(sycl::device & dev, sycl::context & ctx) : d(dev), c(ctx) {}
        _device(sycl::device && dev, sycl::context && ctx) : d(dev), c(ctx) {}
    };

    struct _queue
    {
        sycl::queue q;
        device d;
        _queue(sycl::queue & que) : q(que) {}
        _queue(sycl::queue && que) : q(que) {}
    };

}
}
}

#endif /* SRC_WRAPPERS_ONEAPI_W_ONEAPI_GPU_MANAGEMENT_H_ */
