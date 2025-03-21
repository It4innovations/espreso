
#ifndef SRC_GPU_OPERATIONS_SC_SYMM_HCSX_DDNY_TRIA_H
#define SRC_GPU_OPERATIONS_SC_SYMM_HCSX_DDNY_TRIA_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/permutation_data_new.h"
#include "gpu/operations/trsm_hcsx_ddny_tri.h"
#include "gpu/operations/herk_ddnx_ddny_tri.h"
#include "math/operations/convert_csx_csy_map.h"
#include "gpu/operations/convert_ddnx_ddny.h"
#include "gpu/operations/copy_ddnx_ddnx.h"
#include "gpu/operations/permute_ddnx_ddnx.h"



// TODO: fix naming

namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class sc_symm_hcsx_ddny_tria
{
    // compute schur complement of symmetric matrix using the triangular optimization
public:
    struct config
    {
        typename trsm_hcsx_ddny_tri<T,I>::config cfg_trsm;
        typename herk_ddnx_ddny_tri<T,I>::config cfg_herk;
        char order_X = '_';
        char order_L = '_';
    };
    using Treal = utils::remove_complex_t<T>;
public:
    sc_symm_hcsx_ddny_tria() = default;
    sc_symm_hcsx_ddny_tria(const sc_symm_hcsx_ddny_tria &) = delete;
    sc_symm_hcsx_ddny_tria(sc_symm_hcsx_ddny_tria &&) = delete;
    sc_symm_hcsx_ddny_tria & operator=(const sc_symm_hcsx_ddny_tria &) = delete;
    sc_symm_hcsx_ddny_tria & operator=(sc_symm_hcsx_ddny_tria &&) = delete;
    ~sc_symm_hcsx_ddny_tria() = default;
public:
    void set_config(config cfg_);
    void set_handles(gpu::mgm::queue q_, gpu::spblas::handle spblas_handle_, gpu::dnblas::handle dnblas_handle_);
    void set_coefficients(Treal alpha_);
    void set_h_A11_solver(DirectSparseSolver<T,I> * h_A11_solver_);
    void set_h_A12(MatrixCsxView_new<T,I> * h_A12_);
    // void set_A22(MatrixCsxView_new<T,I> * A22); // assume it is full of zeros
    void set_d_sc(MatrixDenseView_new<T> * d_sc_);
    void setup();
    size_t get_wss_internal();
    size_t get_wss_persistent();
    size_t get_wss_tmp_preprocess();
    size_t get_wss_tmp_perform();
    void set_ws_persistent(void * ws_persistent_);
    void preprocess_submit(void * ws_tmp);
    void perform_submit(void * ws_tmp);
private:
    config cfg;
    gpu::mgm::queue q;
    gpu::spblas::handle handle_spblas;
    gpu::dnblas::handle handle_dnblas;
    DirectSparseSolver<T,I> * h_A11_solver = nullptr;
    MatrixCsxView_new<T,I> * h_A12 = nullptr;
    MatrixDenseView_new<T> * d_sc = nullptr;
    Treal alpha = Treal{1};
    void * ws_persistent = nullptr;
    size_t wss_internal = 0;
    size_t wss_persistent = 0;
    size_t wss_tmp_preprocess = 0;
    size_t wss_tmp_perform = 0;
    bool called_set_config = false;
    bool called_set_handles = false;
    bool called_setup = false;
    bool called_preprocess = false;
private:
    std::unique_ptr<AllocatorArena_new> ator_ws_persistent;
    std::unique_ptr<AllocatorArena_new> ator_ws_tmp_linear;
    std::unique_ptr<AllocatorSinglePointer_new> ator_ws_tmp_overlap;
    size_t wss_tmp_preprocess_linear = 0;
    size_t wss_tmp_preprocess_overlap = 0;
    size_t wss_tmp_perform_linear = 0;
    size_t wss_tmp_perform_overlap = 0;
    size_t A11size = 0;
    bool need_reorder_factor_L2U = false;
    bool need_reorder_factor_U2L = false;
    char solver_factor_uplo = '_';
    MatrixCsxData_new<T,I> h_factor_row_U;
    MatrixCsxData_new<T,I> h_factor_row_L;
    MatrixCsxView_new<T,I> h_L_row;
    MatrixCsxView_new<T,I> h_L_col;
    MatrixCsxView_new<T,I> * h_L = nullptr;
    PermutationData_new<I> h_perm_to_sort_A12_cols;
    PermutationData_new<I> d_perm_to_sort_A12_cols;
    PermutationView_new<I> d_perm_to_sort_back_sc;
    PermutationView_new<I> h_perm_fillreduce;
    MatrixCsxData_new<T,I> h_X_sp;
    MatrixCsxData_new<T,I> d_X_sp;
    MatrixDenseData_new<T> d_X_dn;
    MatrixDenseData_new<T> d_sc_tmp1;
    MatrixDenseData_new<T> d_sc_tmp2;
    std::unique_ptr<convert_dcsx_ddny<T,I>> op_X_sp2dn;
    math::operations::convert_csx_csy_map<T,I> op_L2U;
    math::operations::convert_csx_csy_map<T,I> op_U2L;
    trsm_hcsx_ddny_tri<T,I> op_d_trsm;
    herk_ddnx_ddny_tri<T,I> op_d_herk;
    std::unique_ptr<convert_ddnx_ddny<T>> op_d_sc_trans;
    std::unique_ptr<copy_ddnx_ddnx<T>> op_d_copy_sc_tmp;
    std::unique_ptr<permute_ddnx_ddnx<T,I>> op_d_perm_sc;
    std::unique_ptr<copy_ddnx_ddnx<T>> op_d_copy_sc_final;
};



}
}
}

#endif /* SRC_GPU_OPERATIONS_SC_SYMM_HCSX_DDNY_TRIA_H */
