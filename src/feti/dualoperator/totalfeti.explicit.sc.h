
#ifndef SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_SC_H_
#define SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_SC_H_

#include "dualoperator.h"
#include "math/wrappers/math.sc_solver.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_dnblas.h"

#include "math/math.h"

namespace espreso {

template <typename T, typename I>
class TotalFETIExplicitSc: public DualOperator<T> {
public:
    TotalFETIExplicitSc(FETI<T> &feti, bool apply_on_gpu_ = false);
    ~TotalFETIExplicitSc();

    void info();
    void set(const step::Step &step);
    void update(const step::Step &step);

    // y = F * x
    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y);
    // y = K+(f - Bt * x)
    void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y);

protected:
    void print(const step::Step &step);

    using DualOperator<T>::feti;
    using DualOperator<T>::d;

    void _apply_cpu(const Vector_Dual<T> &x, Vector_Dual<T> &y);
    void _apply_gpu(const Vector_Dual<T> &x, Vector_Dual<T> &y);

    void _apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);

private:
    struct per_domain_stuff
    {
        Matrix_Dense<T,I> F;
        Vector_Dense<T,I> x;
        Vector_Dense<T,I> y;
        Matrix_CSR<T,I> Kreg;
        Matrix_CSR<T,I> Bt;
        Matrix_CSR<T,I> concat_matrix; // Kreg Bt B O
        Matrix_CSR<T,I> null_matrix_A22;
        Matrix_CSR<T,I> null_matrix_A21;
        Vector_Dense<I,I> map_B_transpose;
        Vector_Dense<I,I> map_concat;
        SchurComplementSolver<T,I> sc_solver;
        I n_dofs_interface;
        I n_dofs_domain;
        char F_fill;
        Matrix_Dense<T,I,gpu::mgm::Ad> d_F;
        Vector_Dense<I,I,gpu::mgm::Ad> d_applyg_D2C;
        Vector_Dense<T,I,gpu::mgm::Ad> d_apply_x;
        Vector_Dense<T,I,gpu::mgm::Ad> d_apply_y;
        Vector_Dense<T,I,gpu::mgm::Ad> d_apply_z;
        Vector_Dense<T,I,gpu::mgm::Ad> d_apply_w;
        size_t buffersize_tmp_apply;
    };
    std::vector<per_domain_stuff> domain_data;
    std::vector<Matrix_Dense<T,I>> Fs_allocated;
    std::vector<Matrix_Dense<T,I,gpu::mgm::Ad>> d_Fs_allocated;
    Vector_Dense<T,I,gpu::mgm::Ad> d_applyg_x_cluster;
    Vector_Dense<T,I,gpu::mgm::Ad> d_applyg_y_cluster;
    Vector_Dense<T*,I,gpu::mgm::Ad> d_applyg_xs_pointers;
    Vector_Dense<T*,I,gpu::mgm::Ad> d_applyg_ys_pointers;
    Vector_Dense<I,I,gpu::mgm::Ad> d_applyg_n_dofs_interfaces;
    Vector_Dense<I*,I,gpu::mgm::Ad> d_applyg_D2Cs_pointers;
    gpu::mgm::device device;
    size_t n_queues;
    gpu::mgm::queue main_q;
    std::vector<gpu::mgm::queue> queues;
    std::vector<gpu::dnblas::handle> handles_dense;
    std::vector<void*> d_buffers_dense;
    size_t n_domains;
    bool apply_on_gpu;
    bool gpu_wait_update;
    bool gpu_wait_apply;
    bool gpu_apply_parallel;
    bool is_system_hermitian;
    bool Fs_share_memory;
};

}

#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_SC_H_ */
