
#ifndef SRC_GPU_OPERATIONS_HERK_DDNX_DDNY_H
#define SRC_GPU_OPERATIONS_HERK_DDNX_DDNY_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_dnblas.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T>
class herk_ddnx_ddny
{
    // C = alpha * A * Ah + beta * C
    // C = alpha * Ah * A + beta * C
public:
    enum struct herk_mode { AhA, AAh };
    using Treal = utils::remove_complex_t<T>;
protected:
    herk_ddnx_ddny() = default;
public:
    herk_ddnx_ddny(const herk_ddnx_ddny &) = delete;
    herk_ddnx_ddny(herk_ddnx_ddny &&) = delete;
    herk_ddnx_ddny & operator=(const herk_ddnx_ddny &) = delete;
    herk_ddnx_ddny & operator=(herk_ddnx_ddny &&) = delete;
    virtual ~herk_ddnx_ddny() = 0;
public:
    static std::unique_ptr<herk_ddnx_ddny<T,I>> make();
public:
    void set_handles(gpu::mgm::queue q_, gpu::dnblas::handle dnblas_handle_);
    void set_matrix_A(MatrixDenseView_new<T,I> A_);
    void set_matrix_C(MatrixDenseView_new<T> C_);
    void set_coefficients(Treal alpha_, Treal beta_)
    void set_mode(herk_mode mode_);
    void setup();
    size_t get_wss_tmp_perform();
    void perform_submit(void * ws_tmp);
    static void submit_all(gpu::mgm::queue q, gpu::dnblas::handle handle_dnblas, MatrixDenseView_new<T> A, MatrixDenseView_new<T> C, Treal alpha, Treal beta, herk_mode mode, Allocator_new * ator_gpu);
protected:
    gpu::mgm::queue q;
    gpu::dnblas::handle handle_dnblas;
    MatrixDenseView_new<T> A;
    MatrixDenseView_new<T> C;
    Treal alpha = Treal{1};
    Treal beta = Treal{0};
    herk_mode mode;
    size_t wss_tmp_perform = 0;
    bool called_set_handles = false;
    bool called_set_A = false;
    bool called_set_C = false;
    bool called_set_mode = false;
    bool called_setup = false;
protected:
    virtual void internal_set_matrix_A() {}
    virtual void internal_set_matrix_C() {}
    virtual void internal_setup() {}
    virtual void internal_perform(void * /*ws_tmp*/) {}
};

herk_ddnx_ddny::~herk_ddnx_ddny() = default;



}
}
}



#endif /* SRC_GPU_OPERATIONS_HERK_DDNX_DDNY_H */
