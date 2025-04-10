
#include "math/operations/trsm_csx_dny_tri.h"

#include <memory>

#include "math/operations/pivots_trails_csx.h"
#include "math/operations/auxiliary/tri_partition_trsm.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void trsm_csx_dny_tri<T,I>::set_config(config cfg_)
{
    cfg = cfg_;

    called_set_config = true;
}



template<typename T, typename I>
void trsm_csx_dny_tri<T,I>::set_L(MatrixCsxView_new<T,I> * L_)
{
    L = L_;
}



template<typename T, typename I>
void trsm_csx_dny_tri<T,I>::set_X(MatrixDenseView_new<T> * X_)
{
    X = X_;
}



template<typename T, typename I>
void trsm_csx_dny_tri<T,I>::calc_X_pattern(MatrixCsxView_new<T,I> & X_pattern)
{
    stacktimer::push("trsm_csx_dny_tri::calc_X_pattern");

    if(called_set_pattern) eslog::error("X patern was already set\n");

    X_colpivots.set(X_pattern.ncols, AllocatorCPU_new::get_singleton());
    X_colpivots.alloc();
    pivots_trails_csx<T,I>::do_all(&X_pattern, &X_colpivots, 'C', 'P', 'B');

    X_rowtrails.set(X_pattern.nrows, AllocatorCPU_new::get_singleton());
    X_rowtrails.alloc();
    pivots_trails_csx<T,I>::do_all(&X_pattern, &X_rowtrails, 'R', 'T', 'F');

    stacktimer::pop();

    called_set_pattern = true;
}



template<typename T, typename I>
void trsm_csx_dny_tri<T,I>::preprocess()
{
    stacktimer::push("trsm_csx_dny_tri::preprocess");

    if(!called_set_config) eslog::error("set config was not called\n");
    if(!called_set_pattern) eslog::error("X pattern was not set\n");
    if(called_preprocess) eslog::error("preproces was already called\n");
    if(L == nullptr) eslog::error("matrix L is not set\n");
    if(X == nullptr) eslog::error("matrix X is not set\n");
    if(L->nrows != L->ncols) eslog::error("L has to be square\n");
    if(L->prop.uplo != 'L') eslog::error("matrix L has to have uplo=L\n");
    if(X->nrows != L->nrows) eslog::error("incompatible matrices\n");
    if(X_colpivots.size != X->ncols) eslog::error("wrong colpivots size\n");
    if(X_rowtrails.size != X->nrows) eslog::error("wrong rowtrails size\n");

    for(size_t i = 1; i < X_colpivots.size; i++) {
        if(X_colpivots.vals[i-1] > X_colpivots.vals[i]) {
            eslog::error("X does not have lower triangular structure\n");
        }
    }
    for(size_t i = 1; i < X_rowtrails.size; i++) {
        if(X_rowtrails.vals[i-1] > X_rowtrails.vals[i]) {
            eslog::error("X does not have lower triangular structure\n");
        }
    }

    tri_partition_trsm partitioner;
    char partition_direction = '_';
    if(cfg.strategy == 'F') partition_direction = 'V';
    if(cfg.strategy == 'R') partition_direction = 'H';
    partitioner.set_config(cfg.partition.algorithm, partition_direction, cfg.partition.parameter);
    partitioner.set_system(X->nrows, X->ncols);
    partitioner.set_output_partition(&partition);
    partitioner.setup();
    num_chunks = partitioner.get_num_chunks();
    partition.set(num_chunks + 1, AllocatorCPU_new::get_singleton());
    partition.alloc();
    partitioner.perform();

    if(cfg.strategy == 'R') {
        ops_chunks_splitrhs.resize(num_chunks);
    }
    if(cfg.strategy == 'F') {
        ops_chunks_splifactor.resize(num_chunks);
    }

    splitrhs_first_dense_chunk = num_chunks;
    for(size_t ch = 0; ch < num_chunks; ch++) {
        if(cfg.strategy == 'R') {
            size_t rhs_start = partition.vals[ch];
            size_t rhs_end = partition.vals[ch+1];

            typename trsm_trirhs_chunk_splitrhs<T,I>::config op_config;
            op_config.factor_order_sp = cfg.splitrhs.factor_order_sp;
            if(cfg.splitrhs.spdn_criteria == 'S') {
                op_config.factor_spdn = 'S';
            }
            if(cfg.splitrhs.spdn_criteria == 'D') {
                op_config.factor_spdn = 'D';
            }
            if(cfg.splitrhs.spdn_criteria == 'C') {
                double fraction_treshold = cfg.splitrhs.spdn_param;
                double curr_fraction = ((double)ch + 0.5) / num_chunks;
                op_config.factor_spdn = ((curr_fraction < fraction_treshold) ? 'S' : 'D');
            }
            if(cfg.splitrhs.spdn_criteria == 'Z') {
                double fraction_treshold = cfg.splitrhs.spdn_param;
                double curr_fraction = (double)(rhs_start + rhs_end) / 2.0 / X->ncols;
                op_config.factor_spdn = ((curr_fraction < fraction_treshold) ? 'S' : 'D');
            }
            if(cfg.splitrhs.spdn_criteria == 'T') {
                size_t k_start = X_colpivots.vals[rhs_start];
                size_t k_size = L->nrows - k_start;
                submatrix_csx_csy<T,I> op_sub_L;
                op_sub_L.set_matrix_src(L);
                op_sub_L.set_bounds(k_start, L->nrows, k_start, L->ncols);
                op_sub_L.setup();
                size_t nnz = op_sub_L.get_output_matrix_nnz() + (L->prop.diag == 'U') * k_size;
                size_t nvals = k_size * (k_size + 1) / 2;
                double fraction_treshold = cfg.splitrhs.spdn_param;
                double curr_fraction = (double)nnz / nvals;
                op_config.factor_spdn = ((curr_fraction < fraction_treshold) ? 'S' : 'D');
            }

            if(splitrhs_first_dense_chunk == num_chunks && op_config.factor_spdn == 'D') {
                splitrhs_first_dense_chunk = ch;
            }

            trsm_trirhs_chunk_splitrhs<T,I> & op_chunk = ops_chunks_splitrhs[ch];
            op_chunk.set_config(op_config);
            op_chunk.set_range(rhs_start, rhs_end);
            op_chunk.set_L_sp(L);
            op_chunk.set_L_dn(&L_dn);
            op_chunk.set_X(X);
            op_chunk.set_X_colpivots(&X_colpivots);
        }
        if(cfg.strategy == 'F') {
            size_t k_start = partition.vals[ch];
            size_t k_end = partition.vals[ch+1];

            typename trsm_trirhs_chunk_splitfactor<T,I>::config op_config;
            op_config.trsm_factor_spdn = cfg.splitfactor.trsm_factor_spdn;
            op_config.trsm_factor_order = cfg.splitfactor.trsm_factor_order;
            op_config.gemm_factor_prune = cfg.splitfactor.gemm_factor_prune;
            
            if(cfg.splitfactor.gemm_spdn_criteria == 'S') {
                op_config.gemm_factor_spdn = 'S';
            }
            if(cfg.splitfactor.gemm_spdn_criteria == 'D') {
                op_config.gemm_factor_spdn = 'D';
            }
            if(cfg.splitfactor.gemm_spdn_criteria == 'C') {
                double fraction_treshold = cfg.splitfactor.gemm_spdn_param;
                double curr_fraction = ((double)ch + 0.5) / num_chunks;
                op_config.gemm_factor_spdn = ((curr_fraction < fraction_treshold) ? 'S' : 'D');
            }
            if(cfg.splitfactor.gemm_spdn_criteria == 'Z') {
                double fraction_treshold = cfg.splitfactor.gemm_spdn_param;
                double curr_fraction = (double)(k_start + k_end) / 2.0 / X->nrows;
                op_config.gemm_factor_spdn = ((curr_fraction < fraction_treshold) ? 'S' : 'D');
            }
            if(cfg.splitfactor.gemm_spdn_criteria == 'T') {
                submatrix_csx_csy<T,I> op_sub_L_bot;
                op_sub_L_bot.set_matrix_src(L);
                op_sub_L_bot.set_bounds(k_end, L->nrows, k_start, k_end);
                op_sub_L_bot.setup();
                size_t nnz = op_sub_L_bot.get_output_matrix_nnz();
                if(cfg.splitfactor.gemm_factor_prune == 'N') {
                    size_t nvals = (L->nrows - k_end) * (k_end - k_start);
                    double fraction_treshold = cfg.splitfactor.gemm_spdn_param;
                    double curr_fraction = (double)nnz / nvals;
                    op_config.gemm_factor_spdn = ((curr_fraction < fraction_treshold) ? 'S' : 'D');
                }
                if(cfg.splitfactor.gemm_factor_prune != 'N') {
                    MatrixCsxData_new<T,I> sub_L_bot_test;
                    sub_L_bot_test.set(L->nrows - k_end, k_end - k_start, nnz, 'R', AllocatorCPU_new::get_singleton());
                    sub_L_bot_test.alloc();
                    op_sub_L_bot.set_matrix_dst(&sub_L_bot_test);
                    op_sub_L_bot.perform();
                    pruning_subset_csx<T,I> op_pruning_subset;
                    op_pruning_subset.set_matrix(&sub_L_bot_test);
                    char pm = cfg.splitfactor.gemm_factor_prune;
                    op_pruning_subset.set_pruning_mode(pm == 'R' || pm == 'A', pm == 'C' || pm == 'A');
                    op_pruning_subset.setup();
                    size_t nvals = op_pruning_subset.get_pruned_nrows() * op_pruning_subset.get_pruned_ncols();
                    sub_L_bot_test.clear();
                    double fraction_treshold = cfg.splitfactor.gemm_spdn_param;
                    double curr_fraction = (double)nnz / nvals;
                    op_config.gemm_factor_spdn = ((curr_fraction < fraction_treshold) ? 'S' : 'D');
                }
            }

            if(op_config.gemm_factor_spdn == 'S') {
                op_config.gemm_factor_order = cfg.splitfactor.gemm_factor_order_sp;
            }
            if(op_config.gemm_factor_spdn == 'D') {
                op_config.gemm_factor_order = cfg.splitfactor.gemm_factor_order_dn;
            }

            trsm_trirhs_chunk_splitfactor<T,I> & op_chunk = ops_chunks_splifactor[ch];
            op_chunk.set_config(op_config);
            op_chunk.set_range(k_start, k_end);
            op_chunk.set_L(L);
            op_chunk.set_X(X);
            op_chunk.set_X_rowtrails(&X_rowtrails);
        }
    }

    if(cfg.strategy == 'R') {
        size_t dense_k_start = partition.vals[splitrhs_first_dense_chunk];
        size_t dense_k_size = L->nrows - dense_k_start;
        L_dn.set(dense_k_size, dense_k_size, cfg.splitrhs.factor_order_dn, AllocatorCPU_new::get_singleton());
        L_dn.prop.uplo = L->prop.uplo;
        L_dn.prop.diag = L->prop.diag;

        op_L_sp2dn.set_matrix_src(L);
        op_L_sp2dn.set_matrix_dst(&L_dn);
    }

    for(size_t ch = 0; ch < num_chunks; ch++) {
        if(cfg.strategy == 'R') {
            ops_chunks_splitrhs[ch].preprocess();
        }
        if(cfg.strategy == 'F') {
            ops_chunks_splifactor[ch].preprocess();
        }
    }

    stacktimer::pop();

    called_preprocess = true;
}



template<typename T, typename I>
void trsm_csx_dny_tri<T,I>::perform()
{
    stacktimer::push("trsm_csx_dny_tri::perform");

    if(!called_preprocess) eslog::error("preprocess was not called\n");

    if(cfg.strategy == 'R' && splitrhs_first_dense_chunk != num_chunks) {
        L_dn.alloc();
        op_L_sp2dn.perform_all();
    }

    for(size_t ch = 0; ch < num_chunks; ch++) {
        if(cfg.strategy == 'R') {
            ops_chunks_splitrhs[ch].perform();
        }
        if(cfg.strategy == 'F') {
            ops_chunks_splifactor[ch].perform();
        }
    }

    L_dn.free();

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class trsm_csx_dny_tri<T,I>;

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
