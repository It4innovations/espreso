
#ifndef SRC_GPU_OPERATIONS_CONVERT_DDNX_DDNY_H
#define SRC_GPU_OPERATIONS_CONVERT_DDNX_DDNY_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "gpu/gpu_management.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
class convert_ddnx_ddny
{
    // respect uplo
protected:
    convert_ddnx_ddny() = default;
public:
    convert_ddnx_ddny(const convert_ddnx_ddny &) = delete;
    convert_ddnx_ddny(convert_ddnx_ddny &&) = delete;
    convert_ddnx_ddny & operator=(const convert_ddnx_ddny &) = delete;
    convert_ddnx_ddny & operator=(convert_ddnx_ddny &&) = delete;
    virtual ~convert_ddnx_ddny() = default;
public:
    static std::unique_ptr<convert_ddnx_ddny<T,I>> make();
public:
    void set_handles(gpu::mgm::queue q_, gpu::dnblas::handle handle_dnblas_);
    void set_matrix_src(MatrixDenseView_new<T,I> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T,I> * M_dst_);
    void setup();
    size_t get_wss_tmp_perform();
    void perform_submit(void * ws_tmp);
protected:
    gpu::mgm::queue q;
    gpu::dnblas::handle handle_dnblas;
    MatrixDenseView_new<T> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
    size_t wss_tmp_perform = 0;
    bool called_set_handles = false;
    bool called_setup = false;
protected:
    virtual void internal_setup() {}
    virtual void internal_perform(void * /*ws_tmp*/) {}
};


    
}
}
}

#endif /* SRC_GPU_OPERATIONS_CONVERT_DDNX_DDNY_H */
