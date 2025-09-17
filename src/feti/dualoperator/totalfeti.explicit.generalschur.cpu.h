
#ifndef SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_GENERAL_CPU_H_
#define SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_GENERAL_CPU_H_

#include "dualoperator.h"
#include "math/math.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/operations/schur_csx_dny.h"
#include "feti/dualoperator/dualop_explicit_applicator.h"

namespace espreso {

template <typename T, typename I>
class TotalFETIExplicitGeneralSchurCpu: public DualOperator<T> {
public:
    TotalFETIExplicitGeneralSchurCpu(FETI<T> &feti);
    virtual ~TotalFETIExplicitGeneralSchurCpu();

    void setup() override;
    size_t get_wss_gpu_persistent() override { return total_wss_gpu_persistent; }
    size_t get_wss_gpu_internal() override { return 0; }
    void set_ws_gpu_persistent(void * ws_gpu_persistent_) override { ws_gpu_persistent = ws_gpu_persistent_; }

    void info() override;
    void set(const step::Step &step) override;
    void update(const step::Step &step) override;

    // y = F * x
    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) override;
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y) override;
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y, const std::vector<int> &filter) override;

    // y = K+(f - Bt * x)
    void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y) override;

protected:
    using DualOperator<T>::feti;
    using DualOperator<T>::d;

private:
    using schur_impl_t = typename math::operations::schur_csx_dny<T,I>::implementation_selector;
    struct config
    {
        char order_F = 'R';
        char mainloop_update_split = 'C'; // Combined, Separate
        bool outer_timers = false;
        bool inner_timers = false;
        bool print_config = false;
        char apply_where = 'C';
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
        MatrixDenseView_new<T> F;
        std::unique_ptr<math::operations::schur_csx_dny<T,I>> op_sc;
    };
    config cfg;
    size_t n_domains = 0;
    std::vector<per_domain_stuff> domain_data;
    std::vector<MatrixDenseData_new<T>> Fs_allocated;
    dualop_explicit_applicator<T,I> applicator;
    size_t total_wss_gpu_persistent = 0;
    void * ws_gpu_persistent = nullptr;
};

}

#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_GENERAL_CPU_H_ */
