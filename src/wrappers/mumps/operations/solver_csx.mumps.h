
#ifndef SRC_WRAPPERS_MUMPS_OPERATIONS_SOLVER_CSX_MUMPS_H
#define SRC_WRAPPERS_MUMPS_OPERATIONS_SOLVER_CSX_MUMPS_H

#include "math/operations/solver_csx.h"



namespace espreso {
namespace math {
namespace operations {



#ifdef HAVE_MUMPS


template<typename T, typename I>
struct solver_csx_mumps_data;

template<typename T, typename I>
class solver_csx_mumps : public solver_csx<T,I>
{
public:
    solver_csx_mumps();
    solver_csx_mumps(const solver_csx_mumps &) = delete;
    solver_csx_mumps(solver_csx_mumps &&) = default;
    solver_csx_mumps & operator=(const solver_csx_mumps &) = delete;
    solver_csx_mumps & operator=(solver_csx_mumps &&) = default;
    virtual ~solver_csx_mumps();
protected:
    void internal_factorize_symbolic() override;
    void internal_factorize_numeric() override;
    void internal_get_permutation(PermutationView_new<I> & perm) override;
    void internal_get_factor_L(MatrixCsxView_new<T,I> & L, bool pattern, bool values) override;
    void internal_get_factor_U(MatrixCsxView_new<T,I> & U, bool pattern, bool values) override;
    void internal_solve(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol) override;
    void internal_solve(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol) override;
    void internal_solve(MatrixCsxView_new<T,I> & rhs, MatrixDenseView_new<T> & sol) override;
private:
    std::unique_ptr<solver_csx_mumps_data<T,I>> data;
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
class solver_csx_mumps : public solver_csx<T,I>
{
public:
    solver_csx_mumps() { eslog::error("mumps wrapper is not available\n"); }
    virtual ~solver_csx_mumps() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_MUMPS_OPERATIONS_SOLVER_CSX_MUMPS_H */
