
#ifndef SRC_MATH_OPERATIONS_HERK_DNX_DNY_TRI_H
#define SRC_MATH_OPERATIONS_HERK_DNX_DNY_TRI_H

#include "math/primitives_new/matrix_dense_view_new.h"
#include "math/primitives_new/vector_dense_data_new.h"
#include "math/primitives_new/matrix_csx_view_new.h"
#include "math/wrappers/math.blas.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class herk_dnx_dny_tri
{
public:
    struct config
    {
        char strategy = '_'; // sTairs, sQuares
        char partition_algorithm = '_'; // Uniform, Minimal work
        int partition_parameter = 0; // depends on algorithm
    };
public:
    herk_dnx_dny_tri() = default;
    herk_dnx_dny_tri(const herk_dnx_dny_tri &) = delete;
    herk_dnx_dny_tri(herk_dnx_dny_tri &&) = delete;
    herk_dnx_dny_tri & operator=(const herk_dnx_dny_tri &) = delete;
    herk_dnx_dny_tri & operator=(herk_dnx_dny_tri &&) = delete;
    ~herk_dnx_dny_tri();
public:
    void set_config(config cfg_);
    void set_matrix_A(MatrixDenseView_new<T> * A_);
    void set_matrix_C(MatrixDenseView_new<T> * C_);
    void set_coefficients(T alpha_, T beta_);
    void set_mode(blas::herk_mode mode_);
    void calc_A_pattern(MatrixCsxView_new<T,I> & A_pattern);
    void preprocess();
    void perform();
    void finalize();
private:
    config cfg;
    MatrixDenseView_new<T> * A = nullptr;
    MatrixDenseView_new<T> * C = nullptr;
    size_t num_chunks = 0;
    VectorDenseData_new<size_t> partition;
    VectorDenseData_new<I> A_pivots;
    VectorDenseData_new<I> A_trails;
    size_t n = 0;
    size_t k = 0;
    T alpha = T{1};
    T beta = T{0};
    blas::herk_mode mode;
    bool config_set = false;
    bool mode_set = false;
    bool pattern_set = false;
    bool preproces_called = false;
private:
    void perform_AhA();
    void perform_AAh();
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_HERK_DNX_DNY_TRI_H */
