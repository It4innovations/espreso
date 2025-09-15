
#ifndef SRC_FETI_DUALOPERATOR_TOTALFETI_IMPLICIT_GENERALSPSOVLER_CPU_H
#define SRC_FETI_DUALOPERATOR_TOTALFETI_IMPLICIT_GENERALSPSOVLER_CPU_H

#include "dualoperator.h"
#include "math/operations/solver_csx.h"

namespace espreso {

template <typename T, typename I>
class TotalFETIImplicitGeneralSparseSolverCpu: public DualOperator<T> {
public:
    TotalFETIImplicitGeneralSparseSolverCpu(FETI<T> &feti);
    virtual ~TotalFETIImplicitGeneralSparseSolverCpu();

    void info() override;
    void set(const step::Step &step) override;
    void update(const step::Step &step) override;

    // y = B * K+ * Bt * x
    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) override;
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y) override;
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y, const std::vector<int> &filter) override;

    // y = K+(f - Bt * x)
    void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y) override;

private:
    using solver_impl_t = typename math::operations::solver_csx<T,I>::implementation_selector;
    struct config
    {
        bool parallel_set = true;
        bool parallel_update = true;
        bool parallel_apply = true;
        bool outer_timers = false;
        bool inner_timers = false;
        bool print_config = false;
        solver_impl_t solver_impl = solver_impl_t::autoselect;
    };
    void setup_config(config & cfg, const FETIConfiguration & feti_ecf_config);

protected:
    using DualOperator<T>::feti;
    using DualOperator<T>::d;
private:
    struct per_domain_stuff
    {
        size_t n_dofs_domain;
        size_t n_dofs_interface;
        Matrix_CSR<T,I> Kreg_old;
        MatrixCsxView_new<T,I> Kreg;
        std::unique_ptr<math::operations::solver_csx<T,I>> op_solver;
    };
    config cfg;
    size_t n_domains = 0;
    std::vector<per_domain_stuff> domain_data;
};

}



#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_IMPLICIT_GENERALSPSOVLER_CPU_H */
