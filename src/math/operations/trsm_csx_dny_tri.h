
#ifndef SRC_MATH_OPERATIONS_TRSM_CSX_DNY_TRI_H
#define SRC_MATH_OPERATIONS_TRSM_CSX_DNY_TRI_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/vector_dense_data_new.h"
#include "math/operations/auxiliary/trsm_trirhs_chunk_splitrhs.h"
#include "math/operations/auxiliary/trsm_trirhs_chunk_splitfactor.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class trsm_csx_dny_tri
{
public:
    struct config
    {
        char strategy = '_'; // split Rhs, split Factor
        struct {
            char algorithm = '_'; // Uniform, Minimum work
            int parameter = 0; // meaning depends on algorithm
        } partition;
        struct {
            char factor_order_sp = '_'; // Rowmajor, Colmajor
            char factor_order_dn = '_'; // Rowmajor, Colmajor
            char spdn_criteria = '_'; // Sparse only, Dense only, fraction of Chunks is sparse, fraction of siZe is sparse, densiTy of factor part
            double spdn_param = 0;
        } splitrhs;
        struct {
            char trsm_factor_spdn = '_'; // Sparse, Dense
            char trsm_factor_order = '_'; // Rowmajor, Colmajor
            char gemm_factor_prune = '_'; // No, Rows only, Cols only, All
            char gemm_factor_order_sp = '_'; // Rowmajor, Colmajor
            char gemm_factor_order_dn = '_'; // Rowmajor, Colmajor
            char gemm_spdn_criteria = '_'; // Sparse only, Dense only, fraction of Chunks is sparse, fraction of siZe is sparse, densiTy of factor part
            double gemm_spdn_param = 0;
        } splitfactor;
    };
public:
    trsm_csx_dny_tri() = default;
    trsm_csx_dny_tri(const trsm_csx_dny_tri &) = delete;
    trsm_csx_dny_tri(trsm_csx_dny_tri &&) = delete;
    trsm_csx_dny_tri & operator=(const trsm_csx_dny_tri &) = delete;
    trsm_csx_dny_tri & operator=(trsm_csx_dny_tri &&) = delete;
    ~trsm_csx_dny_tri();
public:
    void set_config(config cfg_);
    void set_L(MatrixCsxView_new<T,I> * L_);
    void set_X(MatrixDenseView_new<T> * X_);
    void calc_X_pattern(MatrixCsxView_new<T,I> & X_pattern);
    void preprocess();
    void perform();
    void finalize();
private:
    config cfg;
    MatrixCsxView_new<T,I> * L = nullptr;
    MatrixDenseView_new<T> * X = nullptr;
    size_t num_chunks = 0;
    VectorDenseData_new<size_t> partition;
    VectorDenseData_new<I> X_colpivots;
    VectorDenseData_new<I> X_rowtrails;
    VectorDenseData_new<trsm_trirhs_chunk_splitrhs<T,I>> ops_chunks_splitrhs;
    VectorDenseData_new<trsm_trirhs_chunk_splitfactor<T,I>> ops_chunks_splifactor;
    bool called_set_config = false;
    bool called_set_pattern = false;
    bool called_preprocess = false;
};



}
}
}



#endif /* SRC_MATH_OPERATIONS_TRSM_CSX_DNY_TRI_H */
