
#include "math/operations/trsm_csx_dny_tri.h"

#include <memory>

#include "config/ecf/operations/trsm_csx_dny_tria.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "math/primitives_new/allocator_new.h"
#include "math/operations/pivots_trails_csx.h"
#include "math/operations/auxiliary/tri_partition_trsm.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
void trsm_csx_dny_tri<T,I>::set_L(MatrixCsxView_new<T,I> * L_)
{
    if(L != nullptr) eslog::error("matrix L is already set\n");

    L = L_;
}



template<typename T, typename I>
void trsm_csx_dny_tri<T,I>::set_X(MatrixDenseView_new<T> * X_)
{
    if(X != nullptr) eslog::error("matrix X is already set\n");

    X = X_;
}



template<typename T, typename I>
void trsm_csx_dny_tri<T,I>::calc_X_pattern(MatrixCsxView_new<T,I> & X_pattern)
{
    stacktimer::push("trsm_csx_dny_tri::calc_X_pattern");

    if(called_set_pattern) eslog::error("X patern was already set\n");
    if(!X_pattern.ator->is_data_accessible_cpu()) eslog::error("matrix X_pattern must be cpu-accessible\n");

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

    if(!called_set_pattern) eslog::error("X pattern was not set\n");
    if(called_preprocess) eslog::error("preproces was already called\n");
    if(L == nullptr) eslog::error("matrix L is not set\n");
    if(X == nullptr) eslog::error("matrix X is not set\n");
    if(!L->ator->is_data_accessible_cpu()) eslog::error("matrix L must be cpu-accessible\n");
    if(!X->ator->is_data_accessible_cpu()) eslog::error("matrix X must be cpu-accessible\n");
    if(L->nrows != L->ncols) eslog::error("L has to be square\n");
    if(L->prop.uplo != 'L') eslog::error("matrix L has to have uplo=L\n");
    if(X->nrows != L->nrows) eslog::error("incompatible matrices\n");
    if(X_colpivots.size != X->ncols) eslog::error("wrong colpivots size\n");
    if(X_rowtrails.size != X->nrows) eslog::error("wrong rowtrails size\n");

    setup_config();

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



template<typename T, typename I>
void trsm_csx_dny_tri<T,I>::setup_config()
{
    using ecf_config = TrsmCsxDnyTriaConfig;
    const ecf_config & ecf = info::ecf->operations.trsm_csx_dny_tria;

    switch(ecf.strategy) {
        case ecf_config::TRSM_TRIA_STRATEGY::AUTO:         cfg.strategy = 'F'; break;
        case ecf_config::TRSM_TRIA_STRATEGY::SPLIT_RHS:    cfg.strategy = 'R'; break;
        case ecf_config::TRSM_TRIA_STRATEGY::SPLIT_FACTOR: cfg.strategy = 'F'; break;
    }

    {
        switch(ecf.partition.algorithm) {
            case ecf_config::PARTITION_ALGORITHM::AUTO:         cfg.partition.algorithm = 'U'; break;
            case ecf_config::PARTITION_ALGORITHM::UNIFORM:      cfg.partition.algorithm = 'U'; break;
            case ecf_config::PARTITION_ALGORITHM::MINIMUM_WORK: cfg.partition.algorithm = 'M'; break;
        }

        char partition_strategy = '_';
        switch(ecf.partition.strategy) {
            case ecf_config::PARTITION_STRATEGY::AUTO:        partition_strategy = 'S'; break;
            case ecf_config::PARTITION_STRATEGY::CHUNK_SIZE:  partition_strategy = 'S'; break;
            case ecf_config::PARTITION_STRATEGY::CHUNK_COUNT: partition_strategy = 'C'; break;
        }

        int chunk_size = ecf.partition.chunk_size;
        if(chunk_size == 0) {
            if(info::mesh->dimension == 2 && cfg.strategy == 'F') chunk_size = 200;
            if(info::mesh->dimension == 2 && cfg.strategy == 'R') chunk_size = 100;
            if(info::mesh->dimension == 3 && cfg.strategy == 'F') chunk_size = 200;
            if(info::mesh->dimension == 3 && cfg.strategy == 'R') chunk_size = 100;
        }

        int chunk_count = utils::replace_if_zero(ecf.partition.chunk_count, 20); // not tested

        if(partition_strategy == 'S') cfg.partition.parameter = -chunk_size;
        if(partition_strategy == 'C') cfg.partition.parameter = chunk_count;
    }

    {
        switch(ecf.split_rhs_config.factor_order_sp) {
            case ecf_config::MATRIX_ORDER::AUTO:      cfg.splitrhs.factor_order_sp = 'R'; break;
            case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.splitrhs.factor_order_sp = 'R'; break;
            case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.splitrhs.factor_order_sp = 'C'; break;
        }

        switch(ecf.split_rhs_config.factor_order_dn) {
            case ecf_config::MATRIX_ORDER::AUTO:      cfg.splitrhs.factor_order_dn = 'C'; break;
            case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.splitrhs.factor_order_dn = 'R'; break;
            case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.splitrhs.factor_order_dn = 'C'; break;
        }

        switch(ecf.split_rhs_config.spdn_criteria) {
            case ecf_config::SPDN_CRITERIA::AUTO:                    cfg.splitrhs.spdn_criteria = 'S'; break;
            case ecf_config::SPDN_CRITERIA::SPARSE_ONLY:             cfg.splitrhs.spdn_criteria = 'S'; break;
            case ecf_config::SPDN_CRITERIA::DENSE_ONLY:              cfg.splitrhs.spdn_criteria = 'D'; break;
            case ecf_config::SPDN_CRITERIA::FRACTION_OF_NUM_CHUNKS:  cfg.splitrhs.spdn_criteria = 'C'; break;
            case ecf_config::SPDN_CRITERIA::FRACTION_OF_FACTOR_SIZE: cfg.splitrhs.spdn_criteria = 'Z'; break;
            case ecf_config::SPDN_CRITERIA::FACTOR_DENSITY:          cfg.splitrhs.spdn_criteria = 'T'; break;
        }
        
        double spdn_param_frac_of_num_chunks = utils::replace_if_zero(ecf.split_rhs_config.spdn_param_frac_of_num_chunks, 0.7); // not tested
        double spdn_param_frac_of_factor_size = utils::replace_if_zero(ecf.split_rhs_config.spdn_param_frac_of_factor_size, 0.7); // not tested
        double spdn_param_factor_density = utils::replace_if_zero(ecf.split_rhs_config.spdn_param_factor_density, 0.1); // not tested

        cfg.splitrhs.spdn_param = 0;
        if(cfg.splitrhs.spdn_criteria == 'C') cfg.splitrhs.spdn_param = spdn_param_frac_of_num_chunks;
        if(cfg.splitrhs.spdn_criteria == 'Z') cfg.splitrhs.spdn_param = spdn_param_frac_of_factor_size;
        if(cfg.splitrhs.spdn_criteria == 'T') cfg.splitrhs.spdn_param = spdn_param_factor_density;
    }

    {
        switch(ecf.split_factor_config.trsm_factor_spdn) {
            case ecf_config::SPDN::AUTO:
                cfg.splitfactor.trsm_factor_spdn = ((info::mesh->dimension == 2) ? 'S' : 'D');
                break;
            case ecf_config::SPDN::SPARSE: cfg.splitfactor.trsm_factor_spdn = 'S'; break;
            case ecf_config::SPDN::DENSE:  cfg.splitfactor.trsm_factor_spdn = 'D'; break;
        }

        switch(ecf.split_factor_config.trsm_factor_order) {
            case ecf_config::MATRIX_ORDER::AUTO:      cfg.splitfactor.trsm_factor_order = 'R'; break;
            case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.splitfactor.trsm_factor_order = 'R'; break;
            case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.splitfactor.trsm_factor_order = 'C'; break;
        }

        switch(ecf.split_factor_config.gemm_factor_pruning) {
            case ecf_config::PRUNING_STRATEGY::AUTO:          cfg.splitfactor.gemm_factor_prune = 'R'; break;
            case ecf_config::PRUNING_STRATEGY::NO_PRUNING:    cfg.splitfactor.gemm_factor_prune = 'N'; break;
            case ecf_config::PRUNING_STRATEGY::ROWS_ONLY:     cfg.splitfactor.gemm_factor_prune = 'R'; break;
            case ecf_config::PRUNING_STRATEGY::COLS_ONLY:     cfg.splitfactor.gemm_factor_prune = 'C'; break;
            case ecf_config::PRUNING_STRATEGY::ROWS_AND_COLS: cfg.splitfactor.gemm_factor_prune = 'A'; break;
        }

        switch(ecf.split_factor_config.gemm_factor_order_sp) {
            case ecf_config::MATRIX_ORDER::AUTO:      cfg.splitfactor.gemm_factor_order_sp = 'R'; break;
            case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.splitfactor.gemm_factor_order_sp = 'R'; break;
            case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.splitfactor.gemm_factor_order_sp = 'C'; break;
        }

        switch(ecf.split_factor_config.gemm_factor_order_dn) {
            case ecf_config::MATRIX_ORDER::AUTO:      cfg.splitfactor.gemm_factor_order_dn = 'R'; break;
            case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.splitfactor.gemm_factor_order_dn = 'R'; break;
            case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.splitfactor.gemm_factor_order_dn = 'C'; break;
        }

        switch(ecf.split_factor_config.gemm_spdn_criteria) {
            case ecf_config::SPDN_CRITERIA::AUTO:
                cfg.splitfactor.gemm_spdn_criteria = ((info::mesh->dimension == 3 && cfg.splitfactor.gemm_factor_prune == 'R') ? 'D' : 'S');
                break;
            case ecf_config::SPDN_CRITERIA::SPARSE_ONLY:             cfg.splitfactor.gemm_spdn_criteria = 'S'; break;
            case ecf_config::SPDN_CRITERIA::DENSE_ONLY:              cfg.splitfactor.gemm_spdn_criteria = 'D'; break;
            case ecf_config::SPDN_CRITERIA::FRACTION_OF_NUM_CHUNKS:  cfg.splitfactor.gemm_spdn_criteria = 'C'; break;
            case ecf_config::SPDN_CRITERIA::FRACTION_OF_FACTOR_SIZE: cfg.splitfactor.gemm_spdn_criteria = 'Z'; break;
            case ecf_config::SPDN_CRITERIA::FACTOR_DENSITY:          cfg.splitfactor.gemm_spdn_criteria = 'T'; break;
        }
        
        double spdn_param_frac_of_num_chunks = utils::replace_if_zero(ecf.split_factor_config.spdn_param_frac_of_num_chunks, 0.7); // not tested
        double spdn_param_frac_of_factor_size = utils::replace_if_zero(ecf.split_factor_config.spdn_param_frac_of_factor_size, 0.7); // not tested
        double spdn_param_factor_density = utils::replace_if_zero(ecf.split_factor_config.spdn_param_factor_density, 0.1); // not tested

        cfg.splitfactor.gemm_spdn_param = 0;
        if(cfg.splitfactor.gemm_spdn_criteria == 'C') cfg.splitfactor.gemm_spdn_param = spdn_param_frac_of_num_chunks;
        if(cfg.splitfactor.gemm_spdn_criteria == 'Z') cfg.splitfactor.gemm_spdn_param = spdn_param_frac_of_factor_size;
        if(cfg.splitfactor.gemm_spdn_criteria == 'T') cfg.splitfactor.gemm_spdn_param = spdn_param_factor_density;
    }
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
