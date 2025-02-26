
#include "math/operations/auxiliary/pivots_trails_csx.h"

#include "math/operations/convert_csx_csy.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void pivots_trails_csx<T,I>::set_mode(char row_col_, char pivots_trails_)
{
    row_col = row_col_;
    pivots_trails = pivots_trails_;
}



template<typename T, typename I>
void pivots_trails_csx<T,I>::set_matrix(MatrixCsxView_new<T,I> * M_)
{
    M = M_;
}



template<typename T, typename I>
void pivots_trails_csx<T,I>::set_output_vector(VectorDenseView_new<size_t> * vec_)
{
    vec = vec_;
}



template<typename T, typename I>
void pivots_trails_csx<T,I>::perform()
{
    if(M == nullptr) eslog::error("matrix is not set\n");
    if(vec == nullptr) eslog::error("vector is not set\n");
    if(row_col == '_' || pivots_trails == '_') eslog::error("mode is not set\n");

    bool need_reorder = (M->order == row_col);
    if(need_reorder) {
        MatrixCsxData_new<T,I> M_reordered;
        M_reordered.set(M->nrows, M->ncols, M->nnz, change_order(M->order), AllocatorCPU_new::get_singleton());
        M_reordered.alloc();

        convert_csx_csy<T,I>::do_all(M, &M_reordered);

        if(pivots_trails == 'P') perform_pivots(M_reordered, *vec);
        if(pivots_trails == 'T') perform_trails(M_reordered, *vec);

        M_reordered.clear();
    }
    else {
        if(pivots_trails == 'P') perform_pivots(*M, *vec);
        if(pivots_trails == 'T') perform_trails(*M, *vec);
    }
}



template<typename T, typename I>
void pivots_trails_csx<T,I>::do_all(MatrixCsxView_new<T,I> * M, VectorDenseView_new<size_t> * vec, char row_col, char pivots_trails)
{
    pivots_trails_csx<T,I> instance;
    instance.set_matrix(M);
    instance.set_output_vector(vec);
    instance.set_mode(row_col, pivots_trails);
    instance.perform();
}



template<typename T, typename I>
void pivots_trails_csx<T,I>::perform_pivots(MatrixCsxView_new<T,I> & A, VectorDenseView_new<size_t> & pivots)
{
    size_t A_primary_size = A.get_size_primary();
    if(A_primary_size != pivots.size) eslog::error("wrong output vector size\n");

    for(size_t i = 0; i < A_primary_size; i++) {
        pivots.vals[i] = A.idxs[A.ptrs[i]];
    }
}



template<typename T, typename I>
void pivots_trails_csx<T,I>::perform_trails(MatrixCsxView_new<T,I> & A, VectorDenseView_new<size_t> & trails)
{
    size_t A_primary_size = A.get_size_primary();
    if(A_primary_size != trails.size) eslog::error("wrong output vector size\n");

    for(size_t i = 0; i < A_primary_size; i++) {
        trails.vals[i] = A.idxs[A.ptrs[i + 1] - 1];
    }
}



#define INSTANTIATE_T_I(T,I) \
template class pivots_trails_csx<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T,int32_t) \
    /* INSTANTIATE_T_I(T,int64_t) */

        #define INSTANTIATE \
        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        /* INSTANTIATE_T(std::complex<double>) */

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}
