
#ifndef SRC_MATH_OPERATIONS_SC_SYMM_CSX_DNY_H
#define SRC_MATH_OPERATIONS_SC_SYMM_CSX_DNY_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/permutation_data_new.h"
#include "math/operations/convert_csx_csy_map.h"
#include "math/operations/trsm_csx_dny_tri.h"
#include "math/operations/herk_dnx_dny_tri.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class sc_symm_csx_dny_tria
{
public:
    struct config
    {
        typename trsm_csx_dny_tri<T,I>::config cfg_trsm;
        typename herk_dnx_dny_tri<T,I>::config cfg_herk;
        char order_X = '_';
        char order_L = '_';
    };
    using Treal = utils::remove_complex_t<T>;
public:
    sc_symm_csx_dny_tria() = default;
    sc_symm_csx_dny_tria(const sc_symm_csx_dny_tria &) = delete;
    sc_symm_csx_dny_tria(sc_symm_csx_dny_tria &&) = default;
    sc_symm_csx_dny_tria & operator=(const sc_symm_csx_dny_tria &) = delete;
    sc_symm_csx_dny_tria & operator=(sc_symm_csx_dny_tria &&) = default;
    ~sc_symm_csx_dny_tria() = default;
public:
    void set_config(config cfg_);
    void set_coefficients(Treal alpha_);
    void set_A11_solver(DirectSparseSolver<T,I> * A11_solver_);
    void set_A12(MatrixCsxView_new<T,I> * A12_);
    // void set_A22(MatrixCsxView_new<T,I> * A22_); // assume it is full of zeros
    void set_sc(MatrixDenseView_new<T> * sc_);
    void preprocess();
    void perform();
private:
    config cfg;
    DirectSparseSolver<T,I> * A11_solver = nullptr;
    MatrixCsxView_new<T,I> * A12 = nullptr;
    MatrixDenseView_new<T> * sc = nullptr;
    Treal alpha = Treal{1};
    bool called_set_config = false;
    bool called_preprocess = false;
private:
    size_t A11size = 0;
    bool need_reorder_factor_L2U = false;
    bool need_reorder_factor_U2L = false;
    char solver_factor_uplo = '_';
    MatrixCsxData_new<T,I> factor_L_row;
    MatrixCsxData_new<T,I> factor_U_row;
    MatrixCsxView_new<T,I> L_row;
    MatrixCsxView_new<T,I> L_col;
    MatrixCsxView_new<T,I> U_row;
    MatrixCsxView_new<T,I> U_col;
    MatrixCsxView_new<T,I> * L_to_use = nullptr;
    MatrixCsxData_new<T,I> X_sp;
    MatrixDenseData_new<T> X_dn;
    MatrixDenseData_new<T> sc_tmp1;
    MatrixDenseData_new<T> sc_tmp2;
    PermutationData_new<I> perm_to_sort_A12_cols;
    PermutationView_new<I> perm_to_sort_back_sc;
    PermutationData_new<I> perm_fillreduce;
    convert_csx_csy_map<T,I> op_L2U;
    convert_csx_csy_map<T,I> op_U2L;
    trsm_csx_dny_tri<T,I> op_trsm;
    herk_dnx_dny_tri<T,I> op_herk;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_SC_SYMM_CSX_DNY_H */
