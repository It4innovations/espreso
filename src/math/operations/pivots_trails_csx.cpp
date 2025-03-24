
#include "math/operations/pivots_trails_csx.h"

#include "math/operations/convert_csx_csy.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void pivots_trails_csx<T,I>::set_mode(char row_col_, char pivots_trails_, char completion_)
{
    row_col = row_col_;
    pivots_trails = pivots_trails_;
    completion = completion_;
}



template<typename T, typename I>
void pivots_trails_csx<T,I>::set_matrix(MatrixCsxView_new<T,I> * M_)
{
    M = M_;
}



template<typename T, typename I>
void pivots_trails_csx<T,I>::set_output_vector(VectorDenseView_new<I> * vec_)
{
    vec = vec_;
}



template<typename T, typename I>
void pivots_trails_csx<T,I>::perform()
{
    stacktimer::push("pivots_trails_csx::perform");

    if(M == nullptr) eslog::error("matrix is not set\n");
    if(vec == nullptr) eslog::error("vector is not set\n");
    if(row_col == '_' || pivots_trails == '_') eslog::error("mode is not set\n");

    bool need_reorder = (M->order != row_col);
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

    I end_val = ((need_reorder) ? M->get_size_primary() : M->get_size_secdary());
    if(completion == 'F') {
        complete_forward(*vec, end_val);
    }
    if(completion == 'B') {
        complete_backward(*vec, end_val);
    }

    stacktimer::pop();
}



template<typename T, typename I>
void pivots_trails_csx<T,I>::do_all(MatrixCsxView_new<T,I> * M, VectorDenseView_new<I> * vec, char row_col, char pivots_trails, char completion)
{
    pivots_trails_csx<T,I> instance;
    instance.set_matrix(M);
    instance.set_output_vector(vec);
    instance.set_mode(row_col, pivots_trails, completion);
    instance.perform();
}



template<typename T, typename I>
void pivots_trails_csx<T,I>::perform_pivots(MatrixCsxView_new<T,I> & A, VectorDenseView_new<I> & pivots)
{
    size_t size_primary = A.get_size_primary();
    size_t size_secdary = A.get_size_secdary();
    if(size_primary != pivots.size) eslog::error("wrong output vector size\n");
    I * ptrs = A.ptrs;
    I * idxs = A.idxs;
    I * piv_vals = pivots.vals;

    for(size_t ip = 0; ip < size_primary; ip++) {
        I start = ptrs[ip];
        I end = ptrs[ip+1];
        if(start == end) {
            piv_vals[ip] = size_secdary;
        }
        else {
            piv_vals[ip] = idxs[start];
        }
    }
}



template<typename T, typename I>
void pivots_trails_csx<T,I>::perform_trails(MatrixCsxView_new<T,I> & A, VectorDenseView_new<I> & trails)
{
    size_t size_primary = A.get_size_primary();
    size_t size_secdary = A.get_size_secdary();
    if(size_primary != trails.size) eslog::error("wrong output vector size\n");
    I * ptrs = A.ptrs;
    I * idxs = A.idxs;
    I * trl_vals = trails.vals;

    for(size_t ip = 0; ip < size_primary; ip++) {
        I start = ptrs[ip];
        I end = ptrs[ip+1];
        if(start == end) {
            trl_vals[ip] = size_secdary;
        }
        else {
            trl_vals[ip] = idxs[end - 1];
        }
    }
}



template<typename T, typename I>
void pivots_trails_csx<T,I>::complete_forward(VectorDenseView_new<I> & vec, I end_val)
{
    size_t size = vec.size;
    I * pivtrl = vec.vals;

    I curr = 0;
    for(size_t ip = 0; ip < size; ip++) {
        if(pivtrl[ip] == end_val) {
            pivtrl[ip] = curr;
        }
        curr = pivtrl[ip];
    }
}



template<typename T, typename I>
void pivots_trails_csx<T,I>::complete_backward(VectorDenseView_new<I> & vec, I end_val)
{
    size_t size = vec.size;
    I * pivtrl = vec.vals;

    I curr = end_val - 1;
    for(size_t ip = size - 1; ip < size; ip--) {
        if(pivtrl[ip] == end_val) {
            pivtrl[ip] = curr;
        }
        curr = pivtrl[ip];
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
        INSTANTIATE_T(std::complex<double>)

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}
