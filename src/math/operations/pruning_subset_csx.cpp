
#include "math/operations/pruning_subset_csx.h"

#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void pruning_subset_csx<T,I>::set_matrix(MatrixCsxView_new<T,I> * M_)
{
    M = M_;
}



template<typename T, typename I>
void pruning_subset_csx<T,I>::set_pruning_mode(bool prune_rows_, bool prune_cols_)
{
    prune_rows = prune_rows_;
    prune_cols = prune_cols_;

    called_set_pruning_mode = true;
}



template<typename T, typename I>
void pruning_subset_csx<T,I>::setup()
{
    stacktimer::push("pruning_subset_csx::setup");

    if(!called_set_pruning_mode) eslog::error("pruning mode is not set\n");
    if(called_setup) eslog::error("setup was already called\n");
    if(M == nullptr) eslog::error("matrix is not set\n");

    size_t size_primary = M->get_size_primary();
    size_t size_secdary = M->get_size_secdary();
    size_t nnz = M->nnz;
    I * ptrs = M->ptrs;
    I * idxs = M->idxs;

    prune_primary = (prune_rows && (M->order == 'R')) || (prune_cols && (M->order == 'C'));
    prune_secdary = (prune_rows && (M->order == 'C')) || (prune_cols && (M->order == 'R'));

    if(prune_primary) {
        pruned_size_primary = 0;
        for(size_t ip = 0; ip < size_primary; ip++) {
            I start = ptrs[ip];
            I end = ptrs[ip+1];
            if(end - start > 0) {
                pruned_size_primary++;
            }
        }
    }
    else {
        pruned_size_primary = size_primary;
    }

    if(prune_secdary) {
        nnz_per_secdary.set(size_secdary, AllocatorCPU_new::get_singleton());
        nnz_per_secdary.alloc();
        std::fill_n(nnz_per_secdary.vals, nnz_per_secdary.size, I{0});
        for(size_t i = 0; i < nnz; i++) {
            I is = idxs[i];
            nnz_per_secdary.vals[is]++;
        }
        pruned_size_secdary = 0;
        for(size_t is = 0; is < size_secdary; is++) {
            if(nnz_per_secdary.vals[is] > 0) {
                pruned_size_secdary++;
            }
        }
    }
    else {
        pruned_size_secdary = size_secdary;
    }

    stacktimer::pop();

    called_setup = true;
}



template<typename T, typename I>
size_t pruning_subset_csx<T,I>::get_pruned_nrows()
{
    if(!called_setup) eslog::error("setup was not called\n");

    if(M->order == 'R') return pruned_size_primary;
    if(M->order == 'C') return pruned_size_secdary;
    eslog::error("invalid order\n");
}



template<typename T, typename I>
size_t pruning_subset_csx<T,I>::get_pruned_ncols()
{
    if(!called_setup) eslog::error("setup was not called\n");

    if(M->order == 'R') return pruned_size_secdary;
    if(M->order == 'C') return pruned_size_primary;
    eslog::error("invalid order\n");
}



template<typename T, typename I>
void pruning_subset_csx<T,I>::set_vector_pruned_rows(VectorDenseView_new<I> * nonempty_rows_)
{
    nonempty_rows = nonempty_rows_;
}



template<typename T, typename I>
void pruning_subset_csx<T,I>::set_vector_pruned_cols(VectorDenseView_new<I> * nonempty_cols_)
{
    nonempty_cols = nonempty_cols_;
}



template<typename T, typename I>
void pruning_subset_csx<T,I>::perform()
{
    stacktimer::push("pruning_subset_csx::perform");

    if(!called_setup) eslog::error("setup was not called\n");
    if(prune_rows && nonempty_rows == nullptr) eslog::error("nonempty rows vector is not set\n");
    if(prune_cols && nonempty_cols == nullptr) eslog::error("nonempty cols vector is not set\n");

    size_t size_primary = M->get_size_primary();
    size_t size_secdary = M->get_size_secdary();
    I * ptrs = M->ptrs;

    VectorDenseView_new<I> * nonempty_primary = nullptr;
    if(M->order == 'R') nonempty_primary = nonempty_rows;
    if(M->order == 'C') nonempty_primary = nonempty_cols;
    if(prune_primary) {
        size_t curr_idx = 0;
        for(size_t ip = 0; ip < size_primary; ip++) {
            I start = ptrs[ip];
            I end = ptrs[ip+1];
            if(end - start > 0) {
                nonempty_primary->vals[curr_idx] = ip;
                curr_idx++;
            }
        }
    }
    else {
        if(nonempty_primary != nullptr) {
            for(size_t ip = 0; ip < size_primary; ip++) {
                nonempty_primary->vals[ip] = ip;
            }
        }
    }

    VectorDenseView_new<I> * nonempty_secdary = nullptr;
    if(M->order == 'R') nonempty_secdary = nonempty_cols;
    if(M->order == 'C') nonempty_secdary = nonempty_rows;
    if(prune_secdary) {
        size_t curr_idx = 0;
        for(size_t is = 0; is < size_secdary; is++) {
            if(nnz_per_secdary.vals[is] > 0) {
                nonempty_secdary->vals[curr_idx] = is;
                curr_idx++;
            }
        }
    }
    else {
        if(nonempty_secdary != nullptr) {
            for(size_t is = 0; is < size_secdary; is++) {
                nonempty_secdary->vals[is] = is;
            }
        }
    }

    stacktimer::pop();
}



template<typename T, typename I>
void pruning_subset_csx<T,I>::finalize()
{
    if(called_setup) {
        nnz_per_secdary.clear();
    }
    called_setup = false;
}



#define INSTANTIATE_T_I(T,I) \
template class pruning_subset_csx<T,I>;

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
