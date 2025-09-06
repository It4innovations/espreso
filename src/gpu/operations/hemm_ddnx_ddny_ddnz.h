
#ifndef SRC_GPU_OPERATIONS_HEMM_DDNX_DDNY_DDNZ_H
#define SRC_GPU_OPERATIONS_HEMM_DDNX_DDNY_DDNZ_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_dnblas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
class hemm_ddnx_ddny_ddnz
{
    // C = alpha * A * B + beta * C
protected:
    hemm_ddnx_ddny_ddnz() = default;
public:
    hemm_ddnx_ddny_ddnz(const hemm_ddnx_ddny_ddnz &) = delete;
    hemm_ddnx_ddny_ddnz(hemm_ddnx_ddny_ddnz &&) = delete;
    hemm_ddnx_ddny_ddnz & operator=(const hemm_ddnx_ddny_ddnz &) = delete;
    hemm_ddnx_ddny_ddnz & operator=(hemm_ddnx_ddny_ddnz &&) = delete;
    virtual ~hemm_ddnx_ddny_ddnz() = default;
public:
    static std::unique_ptr<hemm_ddnx_ddny_ddnz<T>> make();
public:
    void set_handles(gpu::mgm::queue q_, gpu::dnblas::handle dnblas_handle_);
    void set_matrix_A(MatrixDenseView_new<T> * A_);
    void set_matrix_B(MatrixDenseView_new<T> * B_);
    void set_matrix_C(MatrixDenseView_new<T> * C_);
    void set_coefficients(T alpha_, T beta_);
    void setup();
    size_t get_wss_tmp_perform();
    void perform_submit(void * ws_tmp);
protected:
    gpu::mgm::queue q;
    gpu::dnblas::handle handle_dnblas;
    MatrixDenseView_new<T> * A = nullptr;
    MatrixDenseView_new<T> * B = nullptr;
    MatrixDenseView_new<T> * C = nullptr;
    size_t wss_tmp_perform = 0;
    T alpha = T{1};
    T beta = T{0};
    bool called_set_handles = false;
    bool called_setup = false;
protected:
    virtual void internal_setup() {}
    virtual void internal_perform(void * /*ws_tmp*/) {}
};



}
}
}

#endif /* SRC_GPU_OPERATIONS_HEMM_DDNX_DDNY_DDNZ_H */
