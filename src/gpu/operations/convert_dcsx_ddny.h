
#ifndef SRC_GPU_OPERATIONS_CONVERT_CSX_DNY_H
#define SRC_GPU_OPERATIONS_CONVERT_CSX_DNY_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_spblas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class convert_dcsx_ddny
{
    // does not respect uplo and diag
protected:
    convert_dcsx_ddny() = default;
public:
    convert_dcsx_ddny(const convert_dcsx_ddny &) = delete;
    convert_dcsx_ddny(convert_dcsx_ddny &&) = delete;
    convert_dcsx_ddny & operator=(const convert_dcsx_ddny &) = delete;
    convert_dcsx_ddny & operator=(convert_dcsx_ddny &&) = delete;
    virtual ~convert_dcsx_ddny() = default;
public:
    static std::unique_ptr<convert_dcsx_ddny<T,I>> make();
public:
    void set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_);
    void set_matrix_src(MatrixCsxView_new<T,I> * M_src_);
    void set_matrix_dst(MatrixDenseView_new<T> * M_dst_);
    void setup();
    size_t get_wss_internal();
    size_t get_wss_persistent();
    size_t get_wss_tmp_preprocess();
    size_t get_wss_tmp_perform();
    void set_ws_persistent(void * ws_persistent_);
    void preprocess_submit(void * ws_tmp);
    void perform_submit(void * ws_tmp);
    static void submit_all(gpu::mgm::queue q, gpu::spblas::handle handle_spblas, MatrixCsxView_new<T,I> * M_src, MatrixDenseView_new<T> * M_dst, Allocator_new * ator_gpu);
protected:
    gpu::mgm::queue q;
    gpu::spblas::handle spblas_handle;
    MatrixCsxView_new<T,I> * M_src;
    MatrixDenseView_new<T> * M_dst;
    void * ws_persistent = nullptr;
    size_t wss_internal = 0;
    size_t wss_persistent = 0;
    size_t wss_tmp_preprocess = 0;
    size_t wss_tmp_perform = 0;
    bool called_set_handles = false;
    bool called_setup = false;
    bool called_preprocess = false;
protected:
    virtual void internal_setup() {}
    virtual void internal_preprocess(void * /*ws_tmp*/) {}
    virtual void internal_perform(void * /*ws_tmp*/) {}

};



}
}
}



#endif /* SRC_GPU_OPERATIONS_CONVERT_CSX_DNY_H */
