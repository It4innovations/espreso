
#ifndef SRC_FETI_DUALOPERATOR_TOTALFETI_IMPLICIT_GENERALSPSOVLER_CPU_H
#define SRC_FETI_DUALOPERATOR_TOTALFETI_IMPLICIT_GENERALSPSOVLER_CPU_H

#include "dualoperator.h"
#include "math/operations/solver_csx.h"

namespace espreso {

template <typename T, typename I>
class TotalFETIImplicitGeneralSpSolverCpu: public DualOperator<T> {
public:
    TotalFETIImplicitGeneralSpSolverCpu(FETI<T> &feti);
    ~TotalFETIImplicitGeneralSpSolverCpu();

    void info();
    void set(const step::Step &step);
    void update(const step::Step &step);

    // y = B * K+ * Bt * x
    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y);

    // y = K+(f - Bt * x)
    void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y);

protected:
    void print(const step::Step &step);

    using DualOperator<T>::feti;
    using DualOperator<T>::d;

public:
    using solver_is_t = typename math::operations::solver_csx<T,I>::implementation_selector;
    struct config
    {
        bool parallel_set = true;
        bool parallel_update = true;
        bool parallel_apply = true;
        bool outer_timers = false;
        bool inner_timers = false;
        bool print_parameters = false;
        solver_is_t solver_is = solver_is_t::autoselect;
    };
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
