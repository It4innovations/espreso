
#include "math/operations/complete_csx_csx_map.h"

#include "math/primitives_new/allocator_new.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void complete_csx_csx_map<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    if(M_src != nullptr) eslog::error("source matrix is already set\n");

    M_src = M_src_;
}



template<typename T, typename I>
void complete_csx_csx_map<T,I>::set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_)
{
    if(M_dst != nullptr) eslog::error("destination matrix is already set\n");

    M_dst = M_dst_;
}



template<typename T, typename I>
void complete_csx_csx_map<T,I>::set_conj(bool do_conj_)
{
    do_conj = do_conj_;
}



template<typename T, typename I>
size_t complete_csx_csx_map<T,I>::get_dst_nnz()
{
    stacktimer::push("complete_csx_csx_map::get_dst_nnz");

    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(!M_src->ator->is_data_accessible_cpu()) eslog::error("source matrix must be cpu-accessible\n");

    size_t diag_vals = 0;
    I src_size_primary = M_src->get_size_primary();
    for(I ip = 0; ip < src_size_primary; ip++) {
        I start = M_src->ptrs[ip];
        I end = M_src->ptrs[ip+1];
        for(I i = start; i < end; i++) {
            I is = M_src->idxs[i];
            if(ip == is) {
                diag_vals++;
            }
        }
    }
    size_t strictly_lower = M_src->nnz - diag_vals;
    size_t full_nnz = 2 * strictly_lower + diag_vals;

    stacktimer::pop();

    return full_nnz;
}



template<typename T, typename I>
void complete_csx_csx_map<T,I>::perform_pattern()
{
    stacktimer::push("complete_csx_csx_map::perform_pattern");

    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(M_src == M_dst) eslog::error("in-place is not supported\n");
    if(!M_src->ator->is_data_accessible_cpu()) eslog::error("source matrix must be cpu-accessible\n");
    if(!M_dst->ator->is_data_accessible_cpu()) eslog::error("destination matrix must be cpu-accessible\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols) eslog::error("matrix sizes dont match\n");
    if(M_src->nrows != M_src->ncols) eslog::error("matrices must be square\n");
    if(M_src->order != M_dst->order) eslog::error("matrix order does not match\n");
    if(!is_uplo(M_src->prop.uplo)) eslog::error("source matrix must have set uplo\n");
    if(is_uplo(M_dst->prop.uplo)) eslog::error("source matrix must not be uplo\n");

    // using terminology for csr, for csc it is equivalent

    struct rci { I r; I c; I i; };
    std::vector<rci> all_entries;
    all_entries.reserve(M_dst->nnz);
    I src_size_primary = M_src->get_size_primary();
    for(I row = 0; row < src_size_primary; row++) {
        I start = M_src->ptrs[row];
        I end = M_src->ptrs[row+1];
        for(I i = start; i < end; i++) {
            I col = M_src->idxs[i];
            all_entries.push_back({row,col,i});
            if(row != col) {
                all_entries.push_back({col,row,i});
            }
        }
    }
    if(all_entries.size() != M_dst->nnz) eslog::error("wrong dst matrix nnz\n");

    std::sort(all_entries.begin(), all_entries.end(), [](const rci & l, const rci & r){ if(l.r != r.r) return l.r < r.r; return l.c < r.c; });

    map.set(M_dst->nnz, AllocatorCPU_new::get_singleton());
    map.alloc();

    I prev_r = -1;
    for(size_t i = 0; i < all_entries.size(); i++) {
        const rci & entry = all_entries[i];
        while(prev_r != entry.r) {
            prev_r++;
            M_dst->ptrs[prev_r] = i;
        }
        M_dst->idxs[i] = entry.c;
        map.vals[i] = entry.i;
    }
    M_dst->ptrs[M_dst->get_size_primary()] = M_dst->nnz;

    called_perform_pattern = true;

    stacktimer::pop();
}



template<typename T, typename I>
void complete_csx_csx_map<T,I>::perform_values()
{
    stacktimer::push("complete_csx_csx_map::perform_values");

    if(!called_perform_pattern) eslog::error("perform pattern was not called\n");

    for(size_t i = 0; i < M_dst->nnz; i++) {
        M_dst->vals[i] = M_src->vals[map.vals[i]];
    }

    if constexpr(utils::is_complex<T>()) if(do_conj) {
        bool toggle = (M_dst->order == 'R') == (M_src->prop.uplo == 'L');
        I size_dst_primary = M_dst->get_size_primary();
        for(I ip = 0; ip < size_dst_primary; ip++) {
            I start = M_dst->ptrs[ip];
            I end = M_dst->ptrs[ip+1];
            for(I i = start; i < end; i++) {
                I is = M_dst->idxs[i];
                T & val = M_dst->vals[i];
                if((ip < is) == toggle) {
                    val = std::conj(val);
                }
            }
        }
    }

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class complete_csx_csx_map<T,I>;

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
