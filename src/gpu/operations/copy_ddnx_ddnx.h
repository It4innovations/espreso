
#ifndef SRC_GPU_OPERATIONS_COPY_DDNX_DDNX_H
#define SRC_GPU_OPERATIONS_COPY_DDNX_DDNX_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "gpu/gpu_management.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
class copy_ddnx_ddnx
{
    // does not respect uplo, has separate parameter for triangle
protected:
    copy_ddnx_ddnx() = default;
public:
    copy_ddnx_ddnx(const copy_ddnx_ddnx &) = delete;
    copy_ddnx_ddnx(copy_ddnx_ddnx &&) = delete;
    copy_ddnx_ddnx & operator=(const copy_ddnx_ddnx &) = delete;
    copy_ddnx_ddnx & operator=(copy_ddnx_ddnx &&) = delete;
    virtual ~copy_ddnx_ddnx() = default;
public:
    static std::unique_ptr<copy_ddnx_ddnx<T,I>> make();
public:
    void set_handles(gpu::mgm::queue q_);
    void set_matrix_src(MatrixDenseView_new<T,I> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T,I> * M_dst_);
    void set_uplo(char uplo_);
    void setup();
    size_t get_wss_tmp_perform();
    void perform_submit(void * ws_tmp);
protected:
    gpu::mgm::queue q;
    MatrixDenseView_new<T> * M_src = nullptr;
    MatrixDenseView_new<T> * M_dst = nullptr;
    char uplo = '_';
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

#endif /* SRC_GPU_OPERATIONS_COPY_DDNX_DDNX_H */
