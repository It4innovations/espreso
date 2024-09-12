
#ifndef SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_SC_H_
#define SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_SC_H_

#include "dualoperator.h"
#include "math/wrappers/math.sc_solver.h"

#include "math/math.h"

namespace espreso {

template <typename T, typename I>
class TotalFETIExplicitSc: public DualOperator<T> {
public:
    TotalFETIExplicitSc(FETI<T> &feti);
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
    };
    std::vector<per_domain_stuff> domain_data;
    size_t n_domains;
};

}

#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_SC_H_ */
