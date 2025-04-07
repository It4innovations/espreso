
#ifndef SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_SCTRIA_H_
#define SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_SCTRIA_H_

#include "dualoperator.h"
#include "math/wrappers/math.sc_solver.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_dnblas.h"

#include "math/math.h"

#include "math/operations/sc_symm_csx_dny_tria.h"

namespace espreso {

template <typename T, typename I>
class TotalFETIExplicitScTria: public DualOperator<T> {
public:
    TotalFETIExplicitScTria(FETI<T> &feti);
    ~TotalFETIExplicitScTria();

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
        bool parallel_set = true;
        bool parallel_update = true;
        bool parallel_apply = true;
        char mainloop_update_split = '_'; // Combined, Separate
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
        MatrixCsxView_new<T,I> Bt;
        MatrixDenseView_new<T> F;
        Matrix_Dense<T,I> F_old;
        math::operations::sc_symm_csx_dny_tria<T,I> op_sc;
        Vector_Dense<T,I> x;
        Vector_Dense<T,I> y;
        char actual_F_uplo = '_';
    };
    config cfg;
    typename math::operations::sc_symm_csx_dny_tria<T,I>::config op_sc_config;
    size_t n_domains = 0;
    std::vector<per_domain_stuff> domain_data;
    std::vector<MatrixDenseData_new<T>> Fs_allocated;
};

}

#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_SCTRIA_H_ */
