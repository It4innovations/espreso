
#include "math/operations/prune_csx_matx.h"

#include "math/operations/fill_dnx.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
prune_csx_matx<T,I>::~prune_csx_matx()
{
    finalize();
}



template<typename T, typename I>
void prune_csx_matx<T,I>::set_pruning_mode(bool prune_rows_, bool prune_cols_)
{
    prune_rows = prune_rows_;
    prune_cols = prune_cols_;

    called_set_pruning_mode = true;
}



template<typename T, typename I>
void prune_csx_matx<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    M_src = M_src_;
}



template<typename T, typename I>
void prune_csx_matx<T,I>::setup()
{
    if(!called_set_pruning_mode) eslog::error("pruning mode is not set\n");
    if(M_src == nullptr) eslog::error("source matrix is not set\n");

    stacktimer::push("prune_csx_matx::setup");

    if(M_src->order == 'R') {
        prune_primary = prune_rows;
        prune_secdary = prune_cols;
    }
    if(M_src->order == 'C') {
        prune_primary = prune_cols;
        prune_secdary = prune_rows;
    }

    op_pruning_subset.set_matrix(M_src);
    op_pruning_subset.set_pruning_mode(prune_rows, prune_cols);
    op_pruning_subset.setup();
    pruned_nrows = op_pruning_subset.get_pruned_nrows();
    pruned_ncols = op_pruning_subset.get_pruned_ncols();

    stacktimer::pop();

    called_setup = true;
}



template<typename T, typename I>
size_t prune_csx_matx<T,I>::get_dst_matrix_nrows()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return pruned_nrows;
}



template<typename T, typename I>
size_t prune_csx_matx<T,I>::get_dst_matrix_ncols()
{
    if(!called_setup) eslog::error("setup was not called\n");

    return pruned_ncols;
}



template<typename T, typename I>
void prune_csx_matx<T,I>::set_vector_pruned_rows(VectorDenseView_new<I> * pruned_rows_)
{
    pruned_rows_vec = pruned_rows_;
}



template<typename T, typename I>
void prune_csx_matx<T,I>::set_vector_pruned_cols(VectorDenseView_new<I> * pruned_cols_)
{
    pruned_cols_vec = pruned_cols_;
}



template<typename T, typename I>
void prune_csx_matx<T,I>::preprocess()
{
    if(!called_preprocess) eslog::error("setup was not called\n");
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(prune_rows && pruned_rows_vec == nullptr) eslog::error("vector for pruned rows is not set\n");
    if(prune_cols && pruned_cols_vec == nullptr) eslog::error("vector for pruned cols is not set\n");
    if(prune_rows && pruned_rows_vec->size != pruned_nrows) eslog::error("wrong size of pruned rows vector\n");
    if(prune_cols && pruned_cols_vec->size != pruned_ncols) eslog::error("wrong size of pruned cols vector\n");

    stacktimer::push("prune_csx_matx::preprocess");

    if(prune_rows) {
        op_pruning_subset.set_vector_pruned_rows(pruned_rows_vec);
    }
    if(prune_cols) {
        op_pruning_subset.set_vector_pruned_cols(pruned_cols_vec);
    }
    op_pruning_subset.perform();
    op_pruning_subset.finalize();

    if(M_src->order == 'R') {
        pruned_idxs_primary = pruned_rows_vec;
        pruned_idxs_secdary = pruned_cols_vec;
    }
    if(M_src->order == 'C') {
        pruned_idxs_primary = pruned_cols_vec;
        pruned_idxs_secdary = pruned_rows_vec;
    }

    if(prune_secdary) {
        pruned_idxs_secdary_inverse.set(M_src->get_size_secdary(), AllocatorCPU_new::get_singleton());
        pruned_idxs_secdary_inverse.alloc();
        std::fill_n(pruned_idxs_secdary_inverse.vals, pruned_idxs_secdary_inverse.size, ~I{0});
        I * pruned_idxs_secdary_vals = pruned_idxs_secdary->vals;
        for(size_t i = 0; i < pruned_idxs_secdary->size; i++) {
            pruned_idxs_secdary_inverse.vals[pruned_idxs_secdary_vals[i]] = i;
        }
    }

    stacktimer::pop();

    called_preprocess2 = true;
}



template<typename T, typename I>
void prune_csx_matx<T,I>::set_matrix_dst_sp(MatrixCsxView_new<T,I> * M_dst_sp_)
{
    M_dst_sp = M_dst_sp_;
}



template<typename T, typename I>
void prune_csx_matx<T,I>::set_matrix_dst_dn(MatrixDenseView_new<T> * M_dst_dn_)
{
    M_dst_dn = M_dst_dn_;
}



template<typename T, typename I>
void prune_csx_matx<T,I>::perform()
{
    if(!called_preprocess2) eslog::error("preprocess was not called\n");

    stacktimer::push("prune_csx_matx::perform");

    if(M_dst_sp != nullptr) {
        perform_sparse();
    }
    if(M_dst_dn != nullptr) {
        perform_dense();
    }

    stacktimer::pop();
}



template<typename T, typename I>
void prune_csx_matx<T,I>::perform_sparse()
{
    if(M_dst_sp->order != M_src->order) eslog::error("matrix orders dont match\n");
    if(M_dst_sp->nrows != pruned_nrows || M_dst_sp->ncols != pruned_ncols || M_dst_sp->nnz != M_src->nnz) eslog::error("wrong destination matrix size\n");

    if(prune_primary) {
        size_t src_size_primary = M_src->get_size_primary();
        I * src_ptrs = M_src->ptrs;
        I * dst_ptrs = M_dst_sp->ptrs;
        size_t curr_idx = 0;
        for(size_t ip = 0; ip < src_size_primary; ip++) {
            I start = src_ptrs[ip];
            I end = src_ptrs[ip+1];
            if(start != end) {
                dst_ptrs[curr_idx] = start;
                curr_idx++;
            }
        }
        dst_ptrs[curr_idx] = M_src->nnz;
    }
    else {
        std::copy_n(M_src->ptrs, M_src->get_size_primary() + 1, M_dst_sp->ptrs);
    }

    if(prune_secdary) {
        size_t nnz = M_src->nnz;
        I * src_idxs = M_src->idxs;
        I * dst_idxs = M_dst_sp->idxs;
        for(size_t i = 0; i < nnz; i++) {
            dst_idxs[i] = pruned_idxs_secdary_inverse.vals[src_idxs[i]];
        }
    }
    else {
        std::copy_n(M_src->idxs, M_src->nnz, M_dst_sp->idxs);
    }

    std::copy_n(M_src->vals, M_src->nnz, M_dst_sp->vals);
}



template<typename T, typename I>
void prune_csx_matx<T,I>::perform_dense()
{
    if(M_dst_dn->order != M_src->order) eslog::error("matrix orders dont match\n");
    if(M_dst_dn->nrows != pruned_nrows || M_dst_dn->ncols != pruned_ncols) eslog::error("wrong destination matrix size\n");

    size_t src_size_primary = M_src->get_size_primary();
    I * src_ptrs = M_src->ptrs;
    I * src_idxs = M_src->idxs;
    T * src_vals = M_src->vals;
    T * dst_vals = M_dst_dn->vals;
    size_t dst_ld = M_dst_dn->ld;

    fill_dnx<T>::do_all(M_dst_dn, T{0});

    if(prune_primary && !prune_secdary) {
        size_t ipd = 0;
        for(size_t ips = 0; ips < src_size_primary; ips++) {
            I start = src_ptrs[ips];
            I end = src_ptrs[ips+1];
            if(start != end) {
                for(I i = start; i < end; i++) {
                    I is = src_idxs[i];
                    T val = src_vals[i];
                    dst_vals[ipd * dst_ld + is] = val;
                }
                ipd++;
            }
        }
    }
    if(!prune_primary && prune_secdary) {
        for(size_t ip = 0; ip < src_size_primary; ip++) {
            I start = src_ptrs[ip];
            I end = src_ptrs[ip+1];
            if(start != end) {
                for(I i = start; i < end; i++) {
                    I iss = src_idxs[i];
                    I isd = pruned_idxs_secdary_inverse.vals[iss];
                    T val = src_vals[i];
                    dst_vals[ip * dst_ld + isd] = val;
                }
            }
        }
    }
    if(prune_primary && prune_secdary) {
        size_t ipd = 0;
        for(size_t ips = 0; ips < src_size_primary; ips++) {
            I start = src_ptrs[ips];
            I end = src_ptrs[ips+1];
            if(start != end) {
                for(I i = start; i < end; i++) {
                    I iss = src_idxs[i];
                    I isd = pruned_idxs_secdary_inverse.vals[iss];
                    T val = src_vals[i];
                    dst_vals[ipd * dst_ld + isd] = val;
                }
                ipd++;
            }
        }
    }

}



template<typename T, typename I>
void prune_csx_matx<T,I>::finalize()
{
    pruned_idxs_secdary_inverse.clear();
    op_pruning_subset.finalize();

    called_setup = false;
}



#define INSTANTIATE_T_I(T,I) \
template class prune_csx_matx<T,I>;

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
