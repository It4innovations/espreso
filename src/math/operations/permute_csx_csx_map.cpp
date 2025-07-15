
#include "math/operations/permute_csx_csx_map.h"
#include "math/primitives_new/allocator_new.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void permute_csx_csx_map<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    if(M_src != nullptr) eslog::error("matrix M_src is already set\n");

    M_src = M_src_;
}



template<typename T, typename I>
void permute_csx_csx_map<T,I>::set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("matrix M_dst is already set\n");

    M_dst = M_dst_;
}



template<typename T, typename I>
void permute_csx_csx_map<T,I>::set_perm_rows(PermutationView_new<I> * perm_rows_)
{
    if(perm_rows != nullptr) eslog::error("perm_rows is already set\n");

    perm_rows = perm_rows_;
}



template<typename T, typename I>
void permute_csx_csx_map<T,I>::set_perm_cols(PermutationView_new<I> * perm_cols_)
{
    if(perm_cols != nullptr) eslog::error("perm_cols is already set\n");

    perm_cols = perm_cols_;
}



template<typename T, typename I>
void permute_csx_csx_map<T,I>::perform_pattern()
{
    stacktimer::push("permute_csx_csx_map::perform_pattern");

    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!M_src->ator->is_data_accessible_cpu()) eslog::error("source matrix must be cpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_cpu()) eslog::error("destination matrix must be cpu-accessible\n");
    if(M_src->order != M_dst->order) eslog::error("matrix orders dont match\n");
    if(!is_uplo_equal(M_src->prop.uplo, M_dst->prop.uplo)) eslog::error("matrix uplo does not match\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols || M_src->nnz != M_dst->nnz) eslog::error("matrix sizes dont match\n");
    if(perm_rows != nullptr && perm_rows->size != M_src->nrows) eslog::error("wrong row perm size\n");
    if(perm_cols != nullptr && perm_cols->size != M_src->ncols) eslog::error("wrong col perm size\n");

    bool is_uplo = (M_src->prop.uplo == 'L' || M_src->prop.uplo == 'U');
    bool toggle = ((M_dst->order == 'R') == (M_dst->prop.uplo == 'L'));

    struct isd_idx { I isd; I idx; }; // index secondary destination, idx of entry in src matrix
    std::vector<std::vector<isd_idx>> out_data(M_dst->get_size_primary());

    PermutationView_new<I> * perm_primary = (M_src->order == 'R') ? perm_rows : perm_cols;
    PermutationView_new<I> * perm_secdary = (M_src->order == 'R') ? perm_cols : perm_rows;

    for(size_t ips = 0; ips < M_src->get_size_primary(); ips++) {
        I start = M_src->ptrs[ips];
        I end = M_src->ptrs[ips+1];
        I ipd = (perm_primary == nullptr) ? ips : perm_primary->src_to_dst[ips];
        for(I i = start; i < end; i++) {
            I iss = M_src->idxs[i];
            I isd = (perm_secdary == nullptr) ? iss : perm_secdary->src_to_dst[iss];
            if(!is_uplo || ((isd < ipd) == toggle)) out_data[ipd].push_back(isd_idx{isd,i});
            else out_data[isd].push_back(isd_idx{ipd,i});
        }
    }

    if(is_uplo) {
        for(size_t ipd = 0; ipd < M_dst->get_size_primary(); ipd++) {
            std::sort(out_data[ipd].begin(), out_data[ipd].end(), [](const isd_idx & l, const isd_idx & r){return l.isd < r.isd;});
        }
    }

    map_dst_to_src.set(M_src->nnz, AllocatorCPU_new::get_singleton());
    map_dst_to_src.alloc();

    I curr_idx = 0;
    for(size_t ipd = 0; ipd < M_dst->get_size_primary(); ipd++) {
        M_dst->ptrs[ipd] = curr_idx;
        for(size_t j = 0; j < out_data[ipd].size(); j++) {
            M_dst->idxs[curr_idx] = out_data[ipd][j].isd;
            map_dst_to_src.vals[curr_idx] = out_data[ipd][j].idx;
            curr_idx++;
        }
    }
    M_dst->ptrs[M_dst->get_size_primary()] = curr_idx;

    stacktimer::pop();
}



template<typename T, typename I>
void permute_csx_csx_map<T,I>::perform_values()
{
    stacktimer::push("permute_csx_csx_map::perform_values");

    for(size_t i = 0; i < M_src->nnz; i++) {
        M_dst->vals[i] = M_src->vals[map_dst_to_src.vals[i]];
    }

    stacktimer::pop();
}



template<typename T, typename I>
void permute_csx_csx_map<T,I>::perform_all()
{
    perform_pattern();
    perform_values();
}



#define INSTANTIATE_T_I(T,I) \
template class permute_csx_csx_map<T,I>;

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
