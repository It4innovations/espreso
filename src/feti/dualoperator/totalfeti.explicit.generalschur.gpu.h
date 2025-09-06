
#ifndef SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_GENERAL_GPU_H_
#define SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_GENERAL_GPU_H_

#include "dualoperator.h"
#include "math/math.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/operations/submatrix_dnx_dnx_view.h"
#include "gpu/operations/schur_hcsx_ddny.h"
#include "feti/dualoperator/dualop_explicit_applicator.h"

namespace espreso {

template <typename T, typename I>
class TotalFETIExplicitGeneralSchurGpu: public DualOperator<T> {
public:
    TotalFETIExplicitGeneralSchurGpu(FETI<T> &feti);
    ~TotalFETIExplicitGeneralSchurGpu();

    void setup() override;
    size_t get_wss_gpu_persistent() override { return total_wss_gpu_persistent; }
    size_t get_wss_gpu_internal() override { return total_wss_gpu_internal; }
    void set_ws_gpu_persistent(void * ws_gpu_persistent_) override { ws_gpu_persistent = ws_gpu_persistent_; }

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

private:
    using schur_impl_t = typename gpu::operations::schur_hcsx_ddny<T,I>::implementation_selector;
    struct config
    {
        char order_F = 'R';
        bool parallel_set = true;
        bool parallel_update = true;
        bool parallel_apply = true;
        char mainloop_update_split = 'C'; // Combined, Separate
        bool gpu_wait_after_mainloop_update = false;
        bool outer_timers = false;
        bool inner_timers = false;
        bool print_config = false;
        char apply_where = 'G';
        schur_impl_t schur_impl = schur_impl_t::autoselect;
    };
    void setup_config(config & cfg, const FETIConfiguration & feti_ecf_config);

private:
    struct per_domain_stuff
    {
        size_t n_dofs_domain;
        size_t n_dofs_interface;
        Matrix_CSR<T,I> Kreg_old;
        MatrixCsxView_new<T,I> Kreg;
        MatrixCsxView_new<T,I> Bt;
        MatrixDenseView_new<T> d_F;
        Matrix_Dense<T,I,gpu::mgm::Ad> d_F_old;
        std::unique_ptr<gpu::operations::schur_hcsx_ddny<T,I>> op_sc;
        math::operations::submatrix_dnx_dnx_view<T> op_sub_F_from_allocd;
    };
    config cfg;
    size_t n_domains = 0;
    size_t n_queues = 0;
    size_t total_wss_gpu_internal = 0;
    size_t total_wss_gpu_persistent = 0;
    void * ws_gpu_persistent = nullptr;
    std::unique_ptr<AllocatorArena_new> ator_ws_gpu_persistent;
    gpu::mgm::queue & main_q;
    std::vector<gpu::mgm::queue> & queues;
    std::vector<gpu::dnblas::handle> & handles_dense;
    std::vector<gpu::spblas::handle> & handles_sparse;
    std::vector<per_domain_stuff> domain_data;
    std::vector<MatrixDenseData_new<T>> d_Fs_allocated;
    dualop_explicit_applicator<T,I> applicator;
};

}

#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_GENERAL_GPU_H_ */
