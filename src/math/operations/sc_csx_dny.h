
#ifndef SRC_MATH_OPERATIONS_SC_CSX_DNY_H
#define SRC_MATH_OPERATIONS_SC_CSX_DNY_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/permutation_data_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class sc_csx_dny
{
public:
    using Treal = utils::remove_complex_t<T>;
    enum struct implementation_selector
    {
        autoselect,
        triangular,
        mklpardiso,
        spsolver
    };
protected:
    sc_csx_dny() = default;
public:
    sc_csx_dny(const sc_csx_dny &) = delete;
    sc_csx_dny(sc_csx_dny &&) = default;
    sc_csx_dny & operator=(const sc_csx_dny &) = delete;
    sc_csx_dny & operator=(sc_csx_dny &&) = default;
    virtual ~sc_csx_dny() = default;
public:
    static std::unique_ptr<sc_csx_dny<T,I>> make(implementation_selector is = implementation_selector::autoselect);
public:
    void set_coefficients(Treal alpha_);
    void set_matrix(MatrixCsxView_new<T,I> * A11_, MatrixCsxView_new<T,I> * A12_, MatrixCsxView_new<T,I> * A21_, MatrixCsxView_new<T,I> * A22_);
    void set_matrix(MatrixCsxView_new<T,I> * A_, size_t size_sc_);
    void set_sc(MatrixDenseView_new<T> * sc_);
    void set_need_solve_A11(bool need_solve_A11_);
    void preprocess();
    // in some SC implementations, perform has two stages, where first is numerical factorization and then actual SC assembly. And I want to be able to express the differences
    void perform_1();
    void perform_2();
    void perform();
    // some SC implementations already factorize A11, and sometimes the user of this class needs to solve with A11 matrix, so I provide this function
    void solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol);
protected:
    MatrixCsxView_new<T,I> * A11 = nullptr;
    MatrixCsxView_new<T,I> * A12 = nullptr;
    MatrixCsxView_new<T,I> * A21 = nullptr;
    MatrixCsxView_new<T,I> * A22 = nullptr;
    MatrixCsxView_new<T,I> * A = nullptr;
    MatrixDenseView_new<T> * sc = nullptr;
    Treal alpha = Treal{1};
    bool need_solve_A11 = false;
    char called_set_matrix = '_';
    bool called_preprocess = false;
    bool called_perform = false;
protected:
    size_t size_matrix = 0;
    size_t size_sc = 0;
    size_t size_A11 = 0;
    bool is_matrix_hermitian = false;
protected:
    virtual void internal_preprocess() {};
    virtual void internal_perform_1() {};
    virtual void internal_perform_2() {};
    virtual void internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol) {};
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_SC_CSX_DNY_H */
