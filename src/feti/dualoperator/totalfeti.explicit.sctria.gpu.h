
#ifndef SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_SCTRIA_GPU_H_
#define SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_SCTRIA_GPU_H_

#include "dualoperator.h"
#include "math/wrappers/math.sc_solver.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_dnblas.h"

#include "math/math.h"

#include "gpu/operations/sc_symm_hcsx_ddny_tria.h"

namespace espreso {

template <typename T, typename I>
class TotalFETIExplicitScTriaGpu: public DualOperator<T> {
public:
    TotalFETIExplicitScTriaGpu(FETI<T> &feti);
    ~TotalFETIExplicitScTriaGpu();

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

    void _apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);

public:
    struct config
    {
        char order_F = '_';
        bool parallel_set = false;
        bool parallel_update = false;
        bool parallel_apply = false;
        char mainloop_update_split = '_'; // Combined, Separate
        bool gpu_wait_after_mainloop_update = false;
        bool outer_timers = false;
        bool inner_timers = false;
        bool print_parameters = false;
    };
private:
    struct per_domain_stuff
    {
        Matrix_CSR<T,I> Kreg;
        DirectSparseSolver<T> solver_Kreg;
        I n_dofs_domain;
        I n_dofs_interface;
        I n_nz_factor;
        MatrixCsxView_new<T,I> h_Bt;
        MatrixDenseView_new<T> d_F;
        Matrix_Dense<T,I,gpu::mgm::Ad> d_F_old;
        std::unique_ptr<gpu::operations::sc_symm_hcsx_ddny_tria<T,I>> op_sc;
        Vector_Dense<T,I,gpu::mgm::Ad> d_apply_x;
        Vector_Dense<T,I,gpu::mgm::Ad> d_apply_y;
        Vector_Dense<T,I,gpu::mgm::Ad> d_apply_z;
        Vector_Dense<T,I,gpu::mgm::Ad> d_apply_w;
        Vector_Dense<I,I,gpu::mgm::Ad> d_applyg_D2C;
    };
    config cfg;
    typename gpu::operations::sc_symm_hcsx_ddny_tria<T,I>::config op_sc_config;
    gpu::mgm::device device;
    gpu::mgm::queue main_q;
    std::vector<gpu::mgm::queue> queues;
    std::vector<gpu::dnblas::handle> handles_dense;
    std::vector<gpu::spblas::handle> handles_sparse;
    std::vector<per_domain_stuff> domain_data;
    size_t n_domains = 0;
    size_t n_queues = 0;
    size_t total_wss_internal = 0;
    size_t total_wss_persistent = 0;
    void * ws_persistent = nullptr;
    size_t wss_tmp_for_cbmba = 0;
    void * ws_tmp_for_cbmba = nullptr;
    std::unique_ptr<AllocatorArena_new> ator_ws_persistent;
    std::unique_ptr<AllocatorCBMB_new> ator_tmp_cbmba;
    std::vector<MatrixDenseData_new<T>> d_Fs_allocated;
    Vector_Dense<T,I,gpu::mgm::Ad> d_applyg_x_cluster;
    Vector_Dense<T,I,gpu::mgm::Ad> d_applyg_y_cluster;
    Vector_Dense<T*,I,gpu::mgm::Ad> d_applyg_xs_pointers;
    Vector_Dense<T*,I,gpu::mgm::Ad> d_applyg_ys_pointers;
    Vector_Dense<I,I,gpu::mgm::Ad> d_applyg_n_dofs_interfaces;
    Vector_Dense<I*,I,gpu::mgm::Ad> d_applyg_D2Cs_pointers;
};

}

#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_SCTRIA_GPU_H_ */
