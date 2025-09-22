
#ifndef SRC_FETI_DUALOPERATOR_HYBRIDFETI_IMPLICIT_GENERALSPSOVLER_CPU_H
#define SRC_FETI_DUALOPERATOR_HYBRIDFETI_IMPLICIT_GENERALSPSOVLER_CPU_H

#include "dualoperator.h"
#include "math/math.h"
#include "math/operations/solver_csx.h"

namespace espreso {

template <typename T, typename I>
class HybridFETIImplicitGeneralSparseSolverCpu: public DualOperator<T> {
public:
    HybridFETIImplicitGeneralSparseSolverCpu(FETI<T> &feti);
    virtual ~HybridFETIImplicitGeneralSparseSolverCpu();

    void info() override;
    void set(const step::Step &step) override;
    void update(const step::Step &step) override;

    // y = B * K+ * Bt * x
    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) override;
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y) override;
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y, const std::vector<int> &filter) override;
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y, const std::vector<std::vector<int>> &filter) override;

    // y = K+(f - Bt * x)
    void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y) override;
    void BtL(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y) override;

private:
    using solver_impl_t = typename math::operations::solver_csx<T,I>::implementation_selector;
    struct config
    {
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



#endif /* SRC_FETI_DUALOPERATOR_HYBRIDFETI_IMPLICIT_GENERALSPSOVLER_CPU_H */
