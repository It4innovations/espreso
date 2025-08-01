
#ifndef SRC_WRAPPERS_SUITESPARSE_OPERATIONS_SOLVER_CSX_CHOLMOD_H
#define SRC_WRAPPERS_SUITESPARSE_OPERATIONS_SOLVER_CSX_CHOLMOD_H

#include "math/operations/solver_csx.h"



namespace espreso {
namespace math {
namespace operations {



#ifdef HAVE_SUITESPARSE


template<typename T, typename I>
struct solver_csx_cholmod_data;

template<typename T, typename I>
class solver_csx_cholmod : public solver_csx<T,I>
{
public:
    solver_csx_cholmod();
    solver_csx_cholmod(const solver_csx_cholmod &) = delete;
    solver_csx_cholmod(solver_csx_cholmod &&) = default;
    solver_csx_cholmod & operator=(const solver_csx_cholmod &) = delete;
    solver_csx_cholmod & operator=(solver_csx_cholmod &&) = default;
    virtual ~solver_csx_cholmod();
protected:
    void internal_factorize_symbolic() override;
    void internal_factorize_numeric() override;
    void internal_get_permutation(PermutationView_new<I> & perm) override;
    void internal_get_factor_L(MatrixCsxView_new<T,I> & L, bool pattern, bool values) override;
    void internal_get_factor_U(MatrixCsxView_new<T,I> & U, bool pattern, bool values) override;
    void internal_solve(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol) override;
    void internal_solve(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol) override;
private:
    void get_factor_impl(MatrixCsxView_new<T,I> & factor, bool pattern, bool values);
private:
    std::unique_ptr<solver_csx_cholmod_data<T,I>> data;
protected:
    using solver_csx<T,I>::A;
    using solver_csx<T,I>::nnz_L;
    using solver_csx<T,I>::nnz_U;
    using solver_csx<T,I>::need_factors;
    using solver_csx<T,I>::need_solve;
    using solver_csx<T,I>::called_factorize_symbolic;
    using solver_csx<T,I>::called_factorize_numeric;
};



#else



template<typename T, typename I>
class solver_csx_cholmod : public solver_csx<T,I>
{
public:
    solver_csx_cholmod() { eslog::error("suitesparse/cholmod wrapper is not available\n"); }
    virtual ~solver_csx_cholmod() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_SUITESPARSE_OPERATIONS_SOLVER_CSX_CHOLMOD_H */
