
#ifndef SRC_WRAPPERS_ONEAPI_W_ONEAPI_GPU_MANAGEMENT_H_
#define SRC_WRAPPERS_ONEAPI_W_ONEAPI_GPU_MANAGEMENT_H_

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#include <sycl/sycl.hpp>
#pragma clang diagnostic pop
#include <unordered_map>
#include <mutex>



namespace espreso {
namespace gpu {
namespace mgm {

    struct _device
    {
        sycl::device d;
        sycl::context c;
        std::vector<queue> qs;
        std::mutex mtx_alloc;
        size_t mem_allocated = 0;
        std::unordered_map<void*,size_t> alloc_sizes;
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
