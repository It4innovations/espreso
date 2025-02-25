
#include "math/operations/auxiliary/sort_csx_csx_pivtrl.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void sort_csx_csx_pivtrl<T,I>::set_matrix_src(MatrixCsxView_new<T,I> * M_src_)
{
    M_src = M_src_;
}



template<typename T, typename I>
void sort_csx_csx_pivtrl<T,I>::set_matrix_dst(MatrixCsxView_new<T,I> * M_dst_)
{
    M_dst = M_dst_;
}



template<typename T, typename I>
void sort_csx_csx_pivtrl<T,I>::set_mode(char row_col_, char piv_trl_, char asc_dsc_)
{
    row_col = row_col_;
    piv_trl = piv_trl_;
    asc_dsc = asc_dsc_;

    mode_set = true;
}



template<typename T, typename I>
void sort_csx_csx_pivtrl<T,I>::perform()
{
    if(M_src == nullptr) eslog::error("source matrix is not set\n");
    if(M_dst == nullptr) eslog::error("destination matrix is not set\n");
    if(!mode_set) eslog::error("mode is not set\n");
    if(M_src->nrows != M_dst->nrows || M_src->ncols != M_dst->ncols || M_src->nnz != M_dst->nnz) eslog::error("matrix sizes dont match\n");
    if(M_src->order != M_dst->order) eslog::error("matrix orders dont match\n");

    VectorDenseData_new<I> pivots_trails;
    if(row_col == 'R') pivots_trails.set(M_src.nrows, AllocatorCPU_new::get_singleton());
    if(row_col == 'C') pivots_trails.set(M_src.ncols, AllocatorCPU_new::get_singleton());
    pivots_trails.alloc();
    pivots_trails_csx<T,I>::do_all(M_src, pivots_trails, row_col, piv_trl);

    struct rc_pt { I rc; I pt; };
    VectorDenseData_new<rc_pt> rcpts;
    rcpts.set(pivots_trails.size, AllocatorCPU_new::get_singleton());
    rcpts.alloc();
    for(size_t i = 0; i < pivots_trails.size; i++) {
        rcpts.vals[i].rc = i;
        rcpts.vals[i].pt = pivots_trails[i];
    }
    pivots_trails.clear();
    if(asc_dsc == 'A') std::sort(rcpts.vals, rcpts.vals + rcpts.size, [](const rc_pt & l, const rc_pt & r){return l.pt < r.pt;});
    if(asc_dsc == 'D') std::sort(rcpts.vals, rcpts.vals + rcpts.size, [](const rc_pt & l, const rc_pt & r){return l.pt > r.pt;});

    PermutationData_new<I> perm;
    perm.set(pivots_trails.size, AllocatorCPU_new::get_singleton());
    perm.alloc();
    for(size_t i = 0; i < rcpts.size; i++) {
        perm.src_to_dst[i] = rcpts.vals[i].rc;
    }
    rcpts.clear();
    PermutationData_new<I>::invert(perm.src_do_dst, perm.dst_to_src, perm.size);

    if(row_col == 'R') permute_csx_csx<T,I>::do_all(M_src, M_dst, perm, nullptr);
    if(row_col == 'C') permute_csx_csx<T,I>::do_all(M_src, M_dst, nullptr, perm);

    perm.clear();
}



template<typename T, typename I>
void sort_csx_csx_pivtrl<T,I>::do_all(MatrixCsxView_new<T,I> * M_src, MatrixCsxView_new<T,I> * M_dst, char row_col, char piv_trl, char asc_dsc)
{
    sort_csx_csx_pivtrl<T,I> instance;
    instance.set_matrix_src(M_src);
    instance.set_matrix_dst(M_dst);
    instance.set_mode(row_col, piv_trl, asc_dsc);
    instance.perform();
}



}
}
}
