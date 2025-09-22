
#ifndef SRC_FETI_DUALOPERATOR_HYBRIDFETI_EXPLICIT_GENERALSCHUR_GPU_H_
#define SRC_FETI_DUALOPERATOR_HYBRIDFETI_EXPLICIT_GENERALSCHUR_GPU_H_

#include "dualoperator.h"
#include "math/math.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/operations/schur_csx_dny.h"
#include "math/operations/submatrix_dnx_dnx_view.h"
#include "gpu/operations/schur_hcsx_ddny.h"
#include "feti/dualoperator/dualop_explicit_applicator.h"

namespace espreso {

template <typename T, typename I>
class HybridFETIExplicitGeneralSchurGpu: public DualOperator<T> {
public:
    HybridFETIExplicitGeneralSchurGpu(FETI<T> &feti);
    virtual ~HybridFETIExplicitGeneralSchurGpu();

    void setup() override;
    size_t get_wss_gpu_persistent() override { return total_wss_gpu_persistent; }
    size_t get_wss_gpu_internal() override { return total_wss_gpu_internal; }
    void set_ws_gpu_persistent(void * ws_gpu_persistent_) override { ws_gpu_persistent = ws_gpu_persistent_; }

    void info() override;
    void set(const step::Step &step) override;
    void update(const step::Step &step) override;

    // y = (B1 * K+ * B1t) * x + (B1 * K+ * B0) * u + (B1 * R * b)
    // y =              F1 * x + (B1 * K+ * B0) * u + (B1 * R * b)
    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) override;
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y) override;
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y, const std::vector<int> &filter) override;
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y, const std::vector<std::vector<int>> &filter) override;

    // y = K+(f - Bt * x)
    void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y) override;
    void BtL(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y) override;

protected:
    using DualOperator<T>::feti;
    using DualOperator<T>::d;
private:
    using schur_impl_t = typename gpu::operations::schur_hcsx_ddny<T,I>::implementation_selector;
    struct config
    {
        char order_F = 'R';
        char mainloop_update_split = 'C'; // Combined, Separate
        bool gpu_wait_after_mainloop_update = false;
        bool outer_timers = false;
        bool inner_timers = false;
        bool print_config = false;
        char apply_where = 'G';
        bool apply_wait_intermediate = false;
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
        MatrixCsxView_new<T,I> B1t;
        MatrixDenseView_new<T> d_F1;
        std::unique_ptr<gpu::operations::schur_hcsx_ddny<T,I>> op_sc;
        math::operations::submatrix_dnx_dnx_view<T> op_sub_F1_from_allocd;
    };
    config cfg;
    size_t n_domains = 0;
    size_t n_dofs_cluster = 0;
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
    std::vector<MatrixDenseData_new<T>> d_F1s_allocated;
    dualop_explicit_applicator<T,I> F1_applicator;
    VectorDenseData_new<T> apply_x_cluster_2;
    VectorDenseData_new<T> apply_y_cluster_2;
private:
    void _apply_hfeti_stuff(const Vector_Dual<T> &x, Vector_Dual<T> &y);
    void _apply_hfeti_stuff(const Matrix_Dual<T> &x, Matrix_Dual<T> &y);
    void _apply_hfeti_stuff(const Matrix_Dual<T> &x, Matrix_Dual<T> &y, const std::vector<int> &filter);
    void _apply_hfeti_stuff(const Matrix_Dual<T> &x, Matrix_Dual<T> &y, const std::vector<std::vector<int>> &filter);
    void _applyK(std::vector<Vector_Dense<T> > &b, std::vector<Vector_Dense<T> > &x, bool do_Kplus_solve);
    void _computeB0();
    void _computeF0();
    void _computeG0();
    void _computeS0();
    void _compute_beta_mu(std::vector<Vector_Dense<T> > &b);

    std::vector<Vector_Dense<T>> Btx, KplusBtx;
    bool isRegularK;
    std::vector<Matrix_CSR<T>> B0;
    std::vector<std::vector<int>> D2C0;
    std::vector<Matrix_Dense<T>> dKB0, origR1;
    std::vector<int> G0offset;
    std::vector<std::vector<int>> permutationF0;
    Matrix_CSR<T> F0, G0;
    Matrix_Dense<T> S0;
    DirectSparseSolver<T> F0Solver;
    DenseSolver<T> Splus;
    Vector_Dense<T> g, beta, mu;
};

}

#endif /* SRC_FETI_DUALOPERATOR_HYBRIDFETI_EXPLICIT_GENERALSCHUR_GPU_H_ */
