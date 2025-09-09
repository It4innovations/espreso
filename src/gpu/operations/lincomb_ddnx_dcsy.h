
#ifndef SRC_GPU_OPERATIONS_LINCOMB_DDNX_DCSY_H
#define SRC_GPU_OPERATIONS_LINCOMB_DDNX_DCSY_H

#include "math/primitives_new.h"
#include "gpu/gpu_management.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class lincomb_ddnx_dcsy
{
// X = alpha * A + beta * B
// inplace allowed (X=A)
// respect uplo
public:
    using Treal = utils::remove_complex_t<T>;
protected:
    lincomb_ddnx_dcsy() = default;
public:
    lincomb_ddnx_dcsy(const lincomb_ddnx_dcsy &) = delete;
    lincomb_ddnx_dcsy(lincomb_ddnx_dcsy &&) = delete;
    lincomb_ddnx_dcsy & operator=(const lincomb_ddnx_dcsy &) = delete;
    lincomb_ddnx_dcsy & operator=(lincomb_ddnx_dcsy &&) = delete;
    virtual ~lincomb_ddnx_dcsy() = default;
public:
    static std::unique_ptr<lincomb_ddnx_dcsy<T,I>> make();
public:
    void set_handles(gpu::mgm::queue q_);
    void set_matrix_X(MatrixDenseView_new<T> * X_);
    void set_matrix_A(MatrixDenseView_new<T> * A_);
    void set_matrix_B(MatrixCsxView_new<T,I> * B_);
    void set_coefficients(Treal alpha_, Treal beta_);
    void setup();
    size_t get_wss_tmp_perform();
    void perform_submit(void * ws_tmp);
protected:
    gpu::mgm::queue q;
    MatrixDenseView_new<T> * X = nullptr;
    MatrixDenseView_new<T> * A = nullptr;
    MatrixCsxView_new<T,I> * B = nullptr;
    Treal alpha = Treal{1};
    Treal beta = Treal{1};
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



#endif /* SRC_GPU_OPERATIONS_LINCOMB_DDNX_DCSY_H */
