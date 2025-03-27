
#ifndef SRC_GPU_OPERATIONS_TRSM_DCSX_DDNY_DDNY_H
#define SRC_GPU_OPERATIONS_TRSM_DCSX_DDNY_DDNY_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_spblas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class trsm_dcsx_ddny_ddny
{
    // solve A * X = B
    // support in-place, where X=B
protected:
    trsm_dcsx_ddny_ddny() = default;
public:
    trsm_dcsx_ddny_ddny(const trsm_dcsx_ddny_ddny &) = delete;
    trsm_dcsx_ddny_ddny(trsm_dcsx_ddny_ddny &&) = delete;
    trsm_dcsx_ddny_ddny & operator=(const trsm_dcsx_ddny_ddny &) = delete;
    trsm_dcsx_ddny_ddny & operator=(trsm_dcsx_ddny_ddny &&) = delete;
    virtual ~trsm_dcsx_ddny_ddny() = default;
public:
    static std::unique_ptr<trsm_dcsx_ddny_ddny<T,I>> make();
public:
    char get_native_place();
    void set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_);
    void set_matrix_A(MatrixCsxView_new<T,I> * A_);
    void set_matrix_X(MatrixDenseView_new<T> * X_);
    void set_matrix_B(MatrixDenseView_new<T> * B_);
    void setup();
    size_t get_wss_internal();
    size_t get_wss_persistent();
    size_t get_wss_tmp_preprocess();
    size_t get_wss_tmp_perform();
    void set_ws_persistent(void * ws_persistent_);
    void preprocess_submit(void * ws_tmp);
    void perform_submit(void * ws_tmp);
protected:
    gpu::mgm::queue q;
    gpu::spblas::handle handle_spblas;
    MatrixCsxView_new<T,I> * A = nullptr;
    MatrixDenseView_new<T> * X = nullptr;
    MatrixDenseView_new<T> * B = nullptr;
    void * ws_persistent = nullptr;
    size_t wss_internal = 0; // approximate number of bytes allocated internally
    size_t wss_persistent = 0;
    size_t wss_tmp_preprocess = 0;
    size_t wss_tmp_perform = 0;
    char place = '_';
    bool called_set_handles = false;
    bool called_setup = false;
    bool called_preprocess = false;
protected:
    virtual char internal_get_native_place() { return '_'; }
    virtual void internal_setup() {}
    virtual void internal_preprocess(void * /*ws_tmp*/) {}
    virtual void internal_perform(void * /*ws_tmp*/) {}
};



}
}
}

#endif /* SRC_GPU_OPERATIONS_TRSM_DCSX_DDNY_DDNY_H */
