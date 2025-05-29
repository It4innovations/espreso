
#ifndef SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_GENERAL_CPU_H_
#define SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_GENERAL_CPU_H_

#include "dualoperator.h"
#include "math/math.h"
#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/operations/sc_csx_dny.h"
#include "feti/dualoperator/dualop_explicit_applicator.h"

namespace espreso {

template <typename T, typename I>
class TotalFETIExplicitGeneralScCpu: public DualOperator<T> {
public:
    TotalFETIExplicitGeneralScCpu(FETI<T> &feti);
    ~TotalFETIExplicitGeneralScCpu();

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
    using sc_is_t = typename math::operations::sc_csx_dny<T,I>::implementation_selector;
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
        sc_is_t sc_is = sc_is_t::autoselect;
    };
private:
    struct per_domain_stuff
    {
        size_t n_dofs_domain;
        size_t n_dofs_interface;
        Matrix_CSR<T,I> Kreg_old;
        MatrixCsxView_new<T,I> Kreg;
        MatrixCsxView_new<T,I> Bt;
        MatrixDenseView_new<T> F;
        Matrix_Dense<T,I> F_old;
        std::unique_ptr<math::operations::sc_csx_dny<T,I>> op_sc;
        char F_uplo_in_rowmajor = '_';
    };
    config cfg;
    size_t n_domains = 0;
    std::vector<per_domain_stuff> domain_data;
    std::vector<MatrixDenseData_new<T>> Fs_allocated;
    dualop_explicit_applicator<T,I> applicator;
};

}

#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_GENERAL_CPU_H_ */
