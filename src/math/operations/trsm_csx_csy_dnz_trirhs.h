
#ifndef SRC_MATH_OPERATIONS_TRSM_CSX_CSY_DNZ_TRIRHS_H
#define SRC_MATH_OPERATIONS_TRSM_CSX_CSY_DNZ_TRIRHS_H

#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/primitives_new/matrix_dense_view_new.h"



template<typename T, typename I>
class trsm_csx_csy_dnz_trirhs
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
            char gemm_factor_prune = '_'; // Yes, No
            char gemm_factor_order_sp = '_'; // Rowmajor, Colmajor
            char gemm_factor_order_dn = '_'; // Rowmajor, Colmajor
            char gemm_spdn_criteria = '_'; // Sparse only, Dense only, fraction of Chunks is sparse, fraction of siZe is sparse, densiTy of factor part
            double gemm_spdn_param = 0;
        } splitfactor;
    };
public:
    trsm_csx_csy_dnz_trirhs() = default;
    trsm_csx_csy_dnz_trirhs(const trsm_csx_csy_dnz_trirhs &) = delete;
    trsm_csx_csy_dnz_trirhs(trsm_csx_csy_dnz_trirhs &&) = delete;
    trsm_csx_csy_dnz_trirhs & operator=(const trsm_csx_csy_dnz_trirhs &) = delete;
    trsm_csx_csy_dnz_trirhs & operator=(trsm_csx_csy_dnz_trirhs &&) = delete;
    ~trsm_csx_csy_dnz_trirhs();
public:
    void set_config(config cfg_);
    void set_L(MatrixCsxView<T,I> * L_);
    void set_B(MatrixCsxView<T,I> * B_);
    void set_X(MatrixDenseView<T> * X_);
    void preprocess();
    void compute();
    void finalize();
private:
    config cfg;
    MatrixCsxView<T,I> * L;
    MatrixCsxView<T,I> * B;
    MatrixDenseView<T> * X;
    size_t num_chunks = 0;
    VectorDenseView_new<size_t> partition;
    VectorDenseView_new<I> B_colpivots;
    VectorDenseView_new<I> B_rowtrails;
    union trsm_tri_rhs {
        trsm_trirhs_chunk_splitrhs<T,I> splitrhs;
        trsm_trirhs_chunk_splitfactor<T,I> splifactor;
    };
    VectorDenseData_new<trsm_tri_rhs> ops_chunks;
    bool called_set_config = false;
    bool called_preprocess = false;
};



#endif /* SRC_MATH_OPERATIONS_TRSM_CSX_CSY_DNZ_TRIRHS_H */
