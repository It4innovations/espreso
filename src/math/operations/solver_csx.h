
#ifndef SRC_MATH_OPERATIONS_SOLVER_CSX
#define SRC_MATH_OPERATIONS_SOLVER_CSX

#include "math/primitives_new.h"


namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class solver_csx
{
public:
    enum struct implementation_selector
    {
        autoselect,
        mklpardiso,
        suitesparse
    };
protected:
    solver_csx() = default;
public:
    solver_csx(const solver_csx &) = delete;
    solver_csx(solver_csx &&) = default;
    solver_csx & operator=(const solver_csx &) = delete;
    solver_csx & operator=(solver_csx &&) = default;
    virtual ~solver_csx() = default;
public:
    static std::unique_ptr<solver_csx<T,I>> make(implementation_selector is = implementation_selector::autoselect, MatrixCsxView_new<T,I> * matrix = nullptr, bool need_factors = false, bool need_solve = false);
public:
    void set_matrix_A(MatrixCsxView_new<T,I> * A_);
    void set_needs(bool need_factors_, bool need_solve_);
    void factorize_symbolic();
    void factorize_numeric();
    size_t get_factor_nnz_L();
    size_t get_factor_nnz_U();
    size_t get_factor_nnz(char uplo = '_'); // recognizes based on uplo
    void get_permutation(PermutationView_new<I> & perm);
    void get_factor_L(MatrixCsxView_new<T,I> & L, bool pattern = true, bool values = true);
    void get_factor_U(MatrixCsxView_new<T,I> & U, bool pattern = true, bool values = true);
    void get_factor(MatrixCsxView_new<T,I> & factor, bool pattern = true, bool values = true); // recognizes based on uplo
    void solve(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol);
    void solve(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol);
protected:
    virtual void internal_factorize_symbolic() {};
    virtual void internal_factorize_numeric() {};
    virtual void internal_get_permutation(PermutationView_new<I> & perm) {};
    virtual void internal_get_factor_L(MatrixCsxView_new<T,I> & L, bool pattern, bool values) {};
    virtual void internal_get_factor_U(MatrixCsxView_new<T,I> & U, bool pattern, bool values) {};
    virtual void internal_solve(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol) {};
    virtual void internal_solve(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol) {};
protected:
    MatrixCsxView_new<T,I> * A = nullptr;
    size_t nnz_L = 0;
    size_t nnz_U = 0;
    bool need_factors = false;
    bool need_solve = false;
    bool called_factorize_symbolic = false;
    bool called_factorize_numeric = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_SOLVER_CSX */
