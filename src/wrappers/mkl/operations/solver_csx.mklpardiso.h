
#ifndef SRC_WRAPPERS_MKL_OPERATIONS_SOLVER_CSX_MKLPARDISO_H
#define SRC_WRAPPERS_MKL_OPERATIONS_SOLVER_CSX_MKLPARDISO_H

#include "math/operations/solver_csx.h"



namespace espreso {
namespace math {
namespace operations {



#ifdef HAVE_MKL


template<typename T, typename I>
struct solver_csx_mklpardiso_data;

template<typename T, typename I>
class solver_csx_mklpardiso : public solver_csx<T,I>
{
public:
    solver_csx_mklpardiso();
    solver_csx_mklpardiso(const solver_csx_mklpardiso &) = delete;
    solver_csx_mklpardiso(solver_csx_mklpardiso &&) = default;
    solver_csx_mklpardiso & operator=(const solver_csx_mklpardiso &) = delete;
    solver_csx_mklpardiso & operator=(solver_csx_mklpardiso &&) = default;
    virtual ~solver_csx_mklpardiso();
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
    std::unique_ptr<solver_csx_mklpardiso_data<T,I>> data;
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
class solver_csx_mklpardiso : public solver_csx<T,I>
{
public:
    solver_csx_mklpardiso() { eslog::error("suitesparse/cholmod wrapper is not available\n"); }
    virtual ~solver_csx_mklpardiso() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_MKL_OPERATIONS_SOLVER_CSX_MKLPARDISO_H */
