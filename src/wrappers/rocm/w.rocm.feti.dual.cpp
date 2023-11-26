
#include "math/wrappers/math.acc.feti.dual.h"
#include "esinfo/eslog.h"

#include <complex>
#include <vector>
#include <algorithm>

#ifdef HAVE_ROCM
#define MY_ROC

#include "my_common.hpp"
#include "matrices.hpp"
#include "roc_stuff.hpp"
#include "hip_stuff_common.hpp"
#include "cholmod_stuff.hpp"

#include <hip/hip_runtime.h>
#include <rocblas/rocblas.h>
#include <rocsparse/rocsparse.h>

namespace espreso {

const char * tool_str = "roc explicit direct";

template<typename T, typename I>
class my_dual_operator_cluster
{
public:
    using value_type = T;
    using int_type = I;
public:
    my_dual_operator_cluster(const char * magicstring_)
        : magicstring{magicstring_}
    { }
    my_dual_operator_cluster(const my_dual_operator_cluster & other) = delete;
    my_dual_operator_cluster(my_dual_operator_cluster && other) = delete;
    my_dual_operator_cluster & operator=(const my_dual_operator_cluster & other) = delete;
    my_dual_operator_cluster & operator=(my_dual_operator_cluster && other) = delete;
    ~my_dual_operator_cluster();
public:
    void set(const std::vector<MatrixCSR<T,I>> & Kregs, const std::vector<MatrixCSR<T,I>> & Bs);
    void update(const std::vector<MatrixCSR<T,I>> & Kregs);
    void apply(MatrixDense<T,I> & y_cluster, const MatrixDense<T,I> & x_cluster, const std::vector<std::vector<I>> & domain_to_cluster_maps);
    void destroy();
public:
    void get_Fs_copies(std::vector<MatrixDense<T,I>> & h_Fs_out)
    {
        if(stage != 2) throw std::runtime_error("invalid stage when calling get_Fs");

        h_Fs_out.clear();
        h_Fs_out.resize(d_Fs.size());
        for(size_t d = 0; d < d_Fs.size(); d++) h_Fs_out[d].resize(d_Fs[d].nrows, d_Fs[d].ncols, -1, true);
        for(size_t d = 0; d < d_Fs.size(); d++) copy_matrix_submit(h_Fs_out[d], d_Fs[d], hipMemcpyDeviceToHost);
        CHECK(hipDeviceSynchronize());
    }
private:
    static constexpr size_t align_B = 512;
    static constexpr size_t align_elem = align_B / sizeof(T);
    static constexpr int rank_gpu_map[] = {4,5,2,3,6,7,0,1}; // assuming LUMI-G
private:
    const char * magicstring;
    int stage = 0; // 0 uninitialized, 1 after set, 2 after update
    size_t n_domains;
    char * mem_pool_device;
    std::unique_ptr<cbmba_resource> cbmba_res_device;
    std::vector<hipStream_t> streams;
    std::vector<rocsparse_handle> handles_sparse;
    std::vector<rocblas_handle> handles_blas;
    std::vector<MatrixDense_hipd<T,I>> d_Fs;
    std::vector<MatrixCSR_hiph<T,I>> h_Us_sp;
    std::vector<MatrixCSR_hiph<T,I>> h_Ls_sp;
    std::vector<MatrixCSR_hiph<T,I>> h_Bperms_sp;
    std::vector<MatrixCSR_hipd<T,I>> d_Us_sp;
    std::vector<MatrixCSR_hipd<T,I>> d_Ls_sp;
    std::vector<MatrixCSR_hipd<T,I>> d_Bperms_sp;
    std::vector<std::unique_ptr<MatrixDense<T,I,cbmba>>> d_Us_dn;
    std::vector<std::unique_ptr<MatrixDense<T,I,cbmba>>> d_Ls_dn;
    std::vector<std::unique_ptr<MatrixDense<T,I,cbmba>>> d_Xs;
    std::vector<std::unique_ptr<MatrixDense<T,I,cbmba>>> d_Ys;
    std::vector<rocsparse_dnmat_descr> gpu_Fs;
    std::vector<rocsparse_spmat_descr> gpu_Us_sp1;
    std::vector<rocsparse_spmat_descr> gpu_Us_sp2;
    std::vector<rocsparse_spmat_descr> gpu_Ls_sp;
    std::vector<rocsparse_spmat_descr> gpu_Bperms_sp;
    std::vector<rocsparse_dnmat_descr> gpu_Us_dn;
    std::vector<rocsparse_dnmat_descr> gpu_Ls_dn;
    std::vector<rocsparse_dnmat_descr> gpu_Xs;
    std::vector<rocsparse_dnmat_descr> gpu_Xs_tc;
    std::vector<rocsparse_dnmat_descr> gpu_Ys_tc;
    std::vector<std::vector<rocsparse_dnvec_descr>> gpu_Xs_tc_cols;
    std::vector<std::vector<rocsparse_dnvec_descr>> gpu_Ys_tc_cols;
    std::vector<cholmod_common*> cm_commons;
    std::vector<cholmod_factor*> cm_factors;
    std::vector<void*> buffers;
    std::vector<size_t> buffersizes;
    std::vector<MatrixDense_hipd<T,I>> d_apply_xs;
    std::vector<MatrixDense_hipd<T,I>> d_apply_ys;
    std::vector<MatrixDense_hiph<T,I>> h_apply_xs;
    std::vector<MatrixDense_hiph<T,I>> h_apply_ys;
};


template<typename T, typename I>
my_dual_operator_cluster<T,I>::~my_dual_operator_cluster()
{
    if(stage != 0) destroy();
}

template<typename T, typename I>
void my_dual_operator_cluster<T,I>::set(const std::vector<MatrixCSR<T,I>> & Kregs, const std::vector<MatrixCSR<T,I>> & Bs)
{
    if(stage != 0) throw std::runtime_error("invalid stage when calling set");

    n_domains = Kregs.size();
    if(Bs.size() != n_domains) throw std::runtime_error("non-matching number of matrices B");

    char joinsplit    = ms_get_joinsplit   (magicstring);
    char waitcontpar  = ms_get_waitcontpar (magicstring);
    char ordering     = ms_get_ordering    (magicstring);
    char spdnfactor   = ms_get_spdnfactor  (magicstring);
    char transcallman = ms_get_transcallman(magicstring);
    char trsvtrsm     = ms_get_trsvtrsm    (magicstring);
    char gemmsyrk     = ms_get_gemmsyrk    (magicstring);

    const bool need_d_Ls_dn = (spdnfactor == 'D' && gemmsyrk == 'G' && transcallman == 'M');
    const bool need_d_Us_dn = (spdnfactor == 'D');
    const bool need_d_Ls_sp = (need_d_Ls_dn || (spdnfactor == 'S' && transcallman == 'M'));
    const bool need_d_Us_sp = (need_d_Us_dn || (spdnfactor == 'S' && (transcallman == 'C' || gemmsyrk == 'G')));
    const bool need_h_Ls_sp = need_d_Ls_sp;
    const bool need_h_Us_sp = true;
    const bool need_d_Bperms_sp = true;
    const bool need_h_Bperms_sp = true;
    const bool need_d_Xs = true;
    const bool need_d_Ys = (spdnfactor == 'S');

    // set the appropriate gpu closest to the used cores
    {
        int mpi_rank = 0; // todo mpi
        // todo check if mpi_size == n_gpus or mpi_size == 1
        // btw, always using 1 gpu per mpi process
        #pragma omp parallel
        {
            CHECK(hipSetDevice(rank_gpu_map[mpi_rank]));
        }
    }

    streams.resize(n_domains);
    handles_sparse.resize(n_domains);
    handles_blas.resize(n_domains);
    d_Fs.resize(n_domains);
    h_Us_sp.resize(n_domains);
    h_Ls_sp.resize(n_domains);
    h_Bperms_sp.resize(n_domains);
    d_Us_sp.resize(n_domains);
    d_Ls_sp.resize(n_domains);
    d_Bperms_sp.resize(n_domains);
    d_Us_dn.resize(n_domains);
    d_Ls_dn.resize(n_domains);
    d_Xs.resize(n_domains);
    d_Ys.resize(n_domains);
    gpu_Fs.resize(n_domains);
    gpu_Us_sp1.resize(n_domains);
    gpu_Us_sp2.resize(n_domains);
    gpu_Ls_sp.resize(n_domains);
    gpu_Bperms_sp.resize(n_domains);
    gpu_Us_dn.resize(n_domains);
    gpu_Ls_dn.resize(n_domains);
    gpu_Xs.resize(n_domains);
    gpu_Xs_tc.resize(n_domains);
    gpu_Ys_tc.resize(n_domains);
    gpu_Xs_tc_cols.resize(n_domains);
    gpu_Ys_tc_cols.resize(n_domains);
    cm_commons.resize(n_domains, nullptr);
    cm_factors.resize(n_domains, nullptr);
    buffers.resize(n_domains, nullptr);
    buffersizes.resize(n_domains, 0);
    d_apply_xs.resize(n_domains);
    d_apply_ys.resize(n_domains);
    h_apply_xs.resize(n_domains);
    h_apply_ys.resize(n_domains);

    #pragma omp parallel for schedule(static,1) if(joinsplit == 'P')
    for(size_t d = 0; d < n_domains; d++)
    {
        I n_dofs_domain = Bs[d].ncols;
        I n_dofs_interface = Bs[d].nrows;
        I ld_domain = ((n_dofs_domain - 1) / align_elem + 1) * align_elem;
        I ld_interface = ((n_dofs_interface - 1) / align_elem + 1) * align_elem;
        I n_nz_factor;

        // symbolic factorization
        {
            my_cholmod_start_and_init<T,I>(&cm_commons[d], ordering);

            cholmod_sparse cm_Kreg = my_cholmod_sparse_view(Kregs[d], 1);

            cm_factors[d] = my_cholmod_analyze<I>(&cm_Kreg, cm_commons[d]);

            n_nz_factor = my_cholmod_factor_get_nnz<I>(cm_factors[d]);
            if(d == 0) printf("Factor: %llux%llu, %llu non-zeros\n", (long long unsigned)n_dofs_domain, (long long unsigned)n_dofs_domain, (long long unsigned)n_nz_factor);
        }

        // init cuda stream and libraries
        {
            CHECK(hipStreamCreate(&streams[d]));
            CHECK(rocsparse_create_handle(&handles_sparse[d]));
            CHECK(rocblas_create_handle(&handles_blas[d]));
            CHECK(rocsparse_set_stream(handles_sparse[d], streams[d]));
            CHECK(rocblas_set_stream(handles_blas[d], streams[d]));
        }

        // create descriptors
        {
            char * dummyptr = reinterpret_cast<char*>(1);
            CHECK(rocsparse_create_csr_descr(   &gpu_Us_sp1[d], n_dofs_domain,    n_dofs_domain, n_nz_factor, dummyptr++, dummyptr++, dummyptr++, rocsparse_index_type<I>(), rocsparse_index_type<I>(), rocsparse_index_base_zero, rocsparse_data_type<T>()));
            CHECK(rocsparse_create_csr_descr(   &gpu_Us_sp2[d], n_dofs_domain,    n_dofs_domain, n_nz_factor, dummyptr++, dummyptr++, dummyptr++, rocsparse_index_type<I>(), rocsparse_index_type<I>(), rocsparse_index_base_zero, rocsparse_data_type<T>()));
            CHECK(rocsparse_create_csr_descr(    &gpu_Ls_sp[d], n_dofs_domain,    n_dofs_domain, n_nz_factor, dummyptr++, dummyptr++, dummyptr++, rocsparse_index_type<I>(), rocsparse_index_type<I>(), rocsparse_index_base_zero, rocsparse_data_type<T>()));
            CHECK(rocsparse_create_csr_descr(&gpu_Bperms_sp[d], n_dofs_interface, n_dofs_domain, Bs[d].nvals, dummyptr++, dummyptr++, dummyptr++, rocsparse_index_type<I>(), rocsparse_index_type<I>(), rocsparse_index_base_zero, rocsparse_data_type<T>()));
            rocsparse_fill_mode upper = rocsparse_fill_mode_upper;
            rocsparse_fill_mode lower = rocsparse_fill_mode_lower;
            rocsparse_diag_type nonunit = rocsparse_diag_type_non_unit;
            // rocsparse_matrix_type triangular = rocsparse_matrix_type_triangular;
            rocsparse_storage_mode sorted = rocsparse_storage_mode_sorted;
            CHECK(rocsparse_spmat_set_attribute(gpu_Us_sp1[d], rocsparse_spmat_fill_mode, &upper, sizeof(upper)));
            CHECK(rocsparse_spmat_set_attribute(gpu_Us_sp1[d], rocsparse_spmat_diag_type, &nonunit, sizeof(nonunit)));
            // CHECK(rocsparse_spmat_set_attribute(gpu_Us_sp1[d], rocsparse_spmat_matrix_type, &triangular, sizeof(triangular)));
            CHECK(rocsparse_spmat_set_attribute(gpu_Us_sp1[d], rocsparse_spmat_storage_mode, &sorted, sizeof(sorted)));
            CHECK(rocsparse_spmat_set_attribute(gpu_Us_sp2[d], rocsparse_spmat_fill_mode, &upper, sizeof(upper)));
            CHECK(rocsparse_spmat_set_attribute(gpu_Us_sp2[d], rocsparse_spmat_diag_type, &nonunit, sizeof(nonunit)));
            // CHECK(rocsparse_spmat_set_attribute(gpu_Us_sp2[d], rocsparse_spmat_matrix_type, &triangular, sizeof(triangular)));
            CHECK(rocsparse_spmat_set_attribute(gpu_Us_sp2[d], rocsparse_spmat_storage_mode, &sorted, sizeof(sorted)));
            CHECK(rocsparse_spmat_set_attribute(gpu_Ls_sp[d], rocsparse_spmat_fill_mode, &lower, sizeof(lower)));
            CHECK(rocsparse_spmat_set_attribute(gpu_Ls_sp[d], rocsparse_spmat_diag_type, &nonunit, sizeof(nonunit)));
            // CHECK(rocsparse_spmat_set_attribute(gpu_Ls_sp[d], rocsparse_spmat_matrix_type, &triangular, sizeof(triangular)));
            CHECK(rocsparse_spmat_set_attribute(gpu_Ls_sp[d], rocsparse_spmat_storage_mode, &sorted, sizeof(sorted)));
            CHECK(rocsparse_create_dnmat_descr(&gpu_Us_dn[d], n_dofs_domain,    n_dofs_domain,    ld_domain, dummyptr++, rocsparse_data_type<T>(), rocsparse_order_row));
            CHECK(rocsparse_create_dnmat_descr(&gpu_Ls_dn[d], n_dofs_domain,    n_dofs_domain,    ld_domain, dummyptr++, rocsparse_data_type<T>(), rocsparse_order_row));
            CHECK(rocsparse_create_dnmat_descr(&gpu_Xs[d],    n_dofs_interface, n_dofs_domain,    ld_domain, dummyptr++, rocsparse_data_type<T>(), rocsparse_order_row));
            CHECK(rocsparse_create_dnmat_descr(&gpu_Xs_tc[d], n_dofs_domain,    n_dofs_interface, ld_domain, dummyptr++, rocsparse_data_type<T>(), rocsparse_order_column));
            CHECK(rocsparse_create_dnmat_descr(&gpu_Ys_tc[d], n_dofs_domain,    n_dofs_interface, ld_domain, dummyptr++, rocsparse_data_type<T>(), rocsparse_order_column));
            gpu_Xs_tc_cols[d].resize(n_dofs_interface);
            for(I j = 0; j < n_dofs_interface; j++) CHECK(rocsparse_create_dnvec_descr(&gpu_Xs_tc_cols[d][j], n_dofs_domain, dummyptr++, rocsparse_data_type<T>()));
            gpu_Ys_tc_cols[d].resize(n_dofs_interface);
            for(I j = 0; j < n_dofs_interface; j++) CHECK(rocsparse_create_dnvec_descr(&gpu_Ys_tc_cols[d][j], n_dofs_domain, dummyptr++, rocsparse_data_type<T>()));
            CHECK(rocsparse_create_dnmat_descr(&gpu_Fs[d], n_dofs_interface, n_dofs_interface, ld_interface, dummyptr++, rocsparse_data_type<T>(), rocsparse_order_column));
        }

        // buffersize
        {
            size_t max_buffersize = 0;
            size_t buffersize = 0;

            T one = 1;
            T zero = 0;

            if(spdnfactor == 'D') CHECK(rocsparse_sparse_to_dense(handles_sparse[d], gpu_Us_sp1[d],    gpu_Us_dn[d], rocsparse_sparse_to_dense_alg_default, &buffersize, nullptr));
            max_buffersize = std::max(buffersize, max_buffersize);

            if(spdnfactor == 'D' && transcallman == 'M') CHECK(rocsparse_sparse_to_dense(handles_sparse[d], gpu_Ls_sp[d],     gpu_Ls_dn[d], rocsparse_sparse_to_dense_alg_default, &buffersize, nullptr));
            max_buffersize = std::max(buffersize, max_buffersize);

            if(true) CHECK(rocsparse_sparse_to_dense(handles_sparse[d], gpu_Bperms_sp[d], gpu_Xs[d],    rocsparse_sparse_to_dense_alg_default, &buffersize, nullptr));
            max_buffersize = std::max(buffersize, max_buffersize);

            if(trsvtrsm == 'V' && spdnfactor == 'S' && transcallman == 'C') CHECK(rocsparse_spsv(handles_sparse[d], rocsparse_operation_transpose, &one, gpu_Us_sp1[d], gpu_Xs_tc_cols[d][0], gpu_Ys_tc_cols[d][0], rocsparse_data_type<T>(), rocsparse_spsv_alg_default, rocsparse_spsv_stage_buffer_size, &buffersize, nullptr));
            max_buffersize = std::max(buffersize, max_buffersize);

            if(trsvtrsm == 'V' && spdnfactor == 'S' && transcallman == 'M') CHECK(rocsparse_spsv(handles_sparse[d], rocsparse_operation_none,      &one, gpu_Ls_sp[d],  gpu_Xs_tc_cols[d][0], gpu_Ys_tc_cols[d][0], rocsparse_data_type<T>(), rocsparse_spsv_alg_default, rocsparse_spsv_stage_buffer_size, &buffersize, nullptr));
            max_buffersize = std::max(buffersize, max_buffersize);

            if(trsvtrsm == 'M' && spdnfactor == 'S' && transcallman == 'C') CHECK(rocsparse_spsm(handles_sparse[d], rocsparse_operation_transpose, rocsparse_operation_none, &one, gpu_Us_sp1[d], gpu_Xs_tc[d], gpu_Ys_tc[d], rocsparse_data_type<T>(), rocsparse_spsm_alg_default, rocsparse_spsm_stage_buffer_size, &buffersize, nullptr));
            max_buffersize = std::max(buffersize, max_buffersize);

            if(trsvtrsm == 'M' && spdnfactor == 'S' && transcallman == 'M') CHECK(rocsparse_spsm(handles_sparse[d], rocsparse_operation_none,      rocsparse_operation_none, &one, gpu_Ls_sp[d],  gpu_Xs_tc[d], gpu_Ys_tc[d], rocsparse_data_type<T>(), rocsparse_spsm_alg_default, rocsparse_spsm_stage_buffer_size, &buffersize, nullptr));
            max_buffersize = std::max(buffersize, max_buffersize);

            if(gemmsyrk == 'G' && trsvtrsm == 'V' && spdnfactor == 'S') CHECK(rocsparse_spsv(handles_sparse[d], rocsparse_operation_none, &one, gpu_Us_sp2[d], gpu_Ys_tc_cols[d][0], gpu_Xs_tc_cols[d][0], rocsparse_data_type<T>(), rocsparse_spsv_alg_default, rocsparse_spsv_stage_buffer_size, &buffersize, nullptr));
            max_buffersize = std::max(buffersize, max_buffersize);

            if(gemmsyrk == 'G' && trsvtrsm == 'M' && spdnfactor == 'S') CHECK(rocsparse_spsm(handles_sparse[d], rocsparse_operation_none, rocsparse_operation_none, &one, gpu_Us_sp2[d], gpu_Ys_tc[d], gpu_Xs_tc[d], rocsparse_data_type<T>(), rocsparse_spsm_alg_default, rocsparse_spsm_stage_buffer_size, &buffersize, nullptr));
            max_buffersize = std::max(buffersize, max_buffersize);

            if(gemmsyrk == 'G') CHECK(rocsparse_spmm(handles_sparse[d], rocsparse_operation_none, rocsparse_operation_none, &one, gpu_Bperms_sp[d], gpu_Xs_tc[d], &zero, gpu_Fs[d], rocsparse_data_type<T>(), rocsparse_spmm_alg_default, rocsparse_spmm_stage_buffer_size, &buffersize, nullptr));
            max_buffersize = std::max(buffersize, max_buffersize);

            CHECK(rocblas_start_device_memory_size_query(handles_blas[d]));
            T * dummyptrT = reinterpret_cast<T*>(sizeof(T));
            if(trsvtrsm == 'V' && spdnfactor == 'D') CHECK(rocblas_dtrsv(handles_blas[d], rocblas_fill_lower, rocblas_operation_none, rocblas_diagonal_non_unit, n_dofs_domain, dummyptrT, ld_domain, dummyptrT, 1));
            if(trsvtrsm == 'M' && spdnfactor == 'D') CHECK(rocblas_dtrsm(handles_blas[d], rocblas_side_left, rocblas_fill_lower, rocblas_operation_none, rocblas_diagonal_non_unit, n_dofs_domain, n_dofs_interface, &one, dummyptrT, ld_domain, dummyptrT, ld_domain));
            if(trsvtrsm == 'V' && spdnfactor == 'D' && transcallman == 'C') CHECK(rocblas_dtrsv(handles_blas[d], rocblas_fill_lower, rocblas_operation_transpose, rocblas_diagonal_non_unit, n_dofs_domain, dummyptrT, ld_domain, dummyptrT, 1));
            if(trsvtrsm == 'V' && spdnfactor == 'D' && transcallman == 'M') CHECK(rocblas_dtrsv(handles_blas[d], rocblas_fill_upper, rocblas_operation_none,      rocblas_diagonal_non_unit, n_dofs_domain, dummyptrT, ld_domain, dummyptrT, 1));
            if(trsvtrsm == 'M' && spdnfactor == 'D' && transcallman == 'C') CHECK(rocblas_dtrsm(handles_blas[d], rocblas_side_left, rocblas_fill_lower, rocblas_operation_transpose, rocblas_diagonal_non_unit, n_dofs_domain, n_dofs_interface, &one, dummyptrT, ld_domain, dummyptrT, ld_domain));
            if(trsvtrsm == 'M' && spdnfactor == 'D' && transcallman == 'M') CHECK(rocblas_dtrsm(handles_blas[d], rocblas_side_left, rocblas_fill_upper, rocblas_operation_none,      rocblas_diagonal_non_unit, n_dofs_domain, n_dofs_interface, &one, dummyptrT, ld_domain, dummyptrT, ld_domain));
            if(spdnfactor == 'S') CHECK(rocblas_dsyrk(handles_blas[d], rocblas_fill_lower, rocblas_operation_transpose, n_dofs_interface, n_dofs_domain, &one, dummyptrT, ld_domain, &zero, dummyptrT, ld_interface));
            if(spdnfactor == 'D') CHECK(rocblas_dsyrk(handles_blas[d], rocblas_fill_lower, rocblas_operation_transpose, n_dofs_interface, n_dofs_domain, &one, dummyptrT, ld_domain, &zero, dummyptrT, ld_interface));
            CHECK(rocblas_dsymv(handles_blas[d], rocblas_fill_lower, n_dofs_interface, &one, dummyptrT, ld_interface, dummyptrT, 1, &zero, dummyptrT, 1));
            CHECK(rocblas_stop_device_memory_size_query(handles_blas[d], &buffersize));
            max_buffersize = std::max(buffersize, max_buffersize);

            buffersizes[d] = max_buffersize;
        }

        // permanent allocations for the lifetime of the program
        {
            // host pinned memory
            if(need_h_Us_sp)         h_Us_sp[d].resize(n_dofs_domain,    n_dofs_domain, n_nz_factor, true);
            if(need_h_Ls_sp)         h_Ls_sp[d].resize(n_dofs_domain,    n_dofs_domain, n_nz_factor, true);
            if(need_h_Bperms_sp) h_Bperms_sp[d].resize(n_dofs_interface, n_dofs_domain, Bs[d].nvals, true);
            h_apply_xs[d].resize(n_dofs_interface, 1, -1, true);
            h_apply_ys[d].resize(n_dofs_interface, 1, -1, true);

            // device memory
            if(need_d_Us_sp)         d_Us_sp[d].resize(n_dofs_domain,    n_dofs_domain,    n_nz_factor,  true);
            if(need_d_Ls_sp)         d_Ls_sp[d].resize(n_dofs_domain,    n_dofs_domain,    n_nz_factor,  true);
            if(need_d_Bperms_sp) d_Bperms_sp[d].resize(n_dofs_interface, n_dofs_domain,    Bs[d].nvals,  true);
            if(true)                    d_Fs[d].resize(n_dofs_interface, n_dofs_interface, ld_interface, true);
            d_apply_xs[d].resize(n_dofs_interface, 1, -1, true);
            d_apply_ys[d].resize(n_dofs_interface, 1, -1, true);
        }

        // set the pointers inside the descriptors of some matrices
        {
            if(need_d_Us_sp)     CHECK(rocsparse_csr_set_pointers(   gpu_Us_sp1[d],     d_Us_sp[d].rowptrs,     d_Us_sp[d].colidxs,     d_Us_sp[d].vals));
            if(need_d_Us_sp)     CHECK(rocsparse_csr_set_pointers(   gpu_Us_sp2[d],     d_Us_sp[d].rowptrs,     d_Us_sp[d].colidxs,     d_Us_sp[d].vals));
            if(need_d_Ls_sp)     CHECK(rocsparse_csr_set_pointers(    gpu_Ls_sp[d],     d_Ls_sp[d].rowptrs,     d_Ls_sp[d].colidxs,     d_Ls_sp[d].vals));
            if(need_d_Bperms_sp) CHECK(rocsparse_csr_set_pointers(gpu_Bperms_sp[d], d_Bperms_sp[d].rowptrs, d_Bperms_sp[d].colidxs, d_Bperms_sp[d].vals));
            rocsparse_dnmat_set_values(gpu_Fs[d], d_Fs[d].vals);
        }

        // prepare matrices on host
        {
            Permutation<I> perm(n_dofs_domain, true);
            my_cholmod_extract_factor_perm(cm_factors[d], perm, cm_commons[d]);
            permute_matrix_cols(h_Bperms_sp[d], Bs[d], perm);
            perm.deallocate();
            if(need_h_Ls_sp) matrix_transpose(h_Ls_sp[d], h_Us_sp[d]);
        }

        // copy some matrices to device
        {
            if(need_d_Bperms_sp) copy_matrix_submit(d_Bperms_sp[d], h_Bperms_sp[d], hipMemcpyHostToDevice, streams[d], true, true, true);
            if(waitcontpar == 'W') CHECK(hipStreamSynchronize(streams[d]));
        }
    }

    // clean up the mess from buggy openmp in clang
    run_dummy_parallel_region();

    // memory pool alloc
    {
        size_t pool_size_device;
        hip_malloc_max_memory(&mem_pool_device, &pool_size_device);
        size_t pool_size_device_aligned = (pool_size_device / align_B) * align_B;
        cbmba_res_device = std::make_unique<cbmba_resource>(mem_pool_device, pool_size_device_aligned);
    }

    CHECK(hipDeviceSynchronize());

    stage = 1;
}

template<typename T, typename I>
void my_dual_operator_cluster<T,I>::update(const std::vector<MatrixCSR<T,I>> & Kregs)
{
    if(stage != 1 && stage != 2) throw std::runtime_error("invalid stage when calling update");

    if(Kregs.size() != n_domains) throw std::runtime_error("non-matching number of matrices Kreg between set and update");

    char joinsplit    = ms_get_joinsplit   (magicstring);
    char waitcontpar  = ms_get_waitcontpar (magicstring);
    char ordering     = ms_get_ordering    (magicstring);
    char spdnfactor   = ms_get_spdnfactor  (magicstring);
    char transcallman = ms_get_transcallman(magicstring);
    char trsvtrsm     = ms_get_trsvtrsm    (magicstring);
    char gemmsyrk     = ms_get_gemmsyrk    (magicstring);

    const bool need_d_Ls_dn = (spdnfactor == 'D' && gemmsyrk == 'G' && transcallman == 'M');
    const bool need_d_Us_dn = (spdnfactor == 'D');
    const bool need_d_Ls_sp = (need_d_Ls_dn || (spdnfactor == 'S' && transcallman == 'M'));
    const bool need_d_Us_sp = (need_d_Us_dn || (spdnfactor == 'S' && (transcallman == 'C' || gemmsyrk == 'G')));
    const bool need_h_Ls_sp = need_d_Ls_sp;
    const bool need_h_Us_sp = true;
    const bool need_d_Bperms_sp = true;
    const bool need_h_Bperms_sp = true;
    const bool need_d_Xs = true;
    const bool need_d_Ys = (spdnfactor == 'S');

    #pragma omp parallel for schedule(static,1) if(waitcontpar == 'P')
    for(size_t d = 0; d < n_domains; d++)
    {
        size_t n_dofs_domain = h_Bperms_sp[d].ncols;
        size_t n_dofs_interface = h_Bperms_sp[d].nrows;
        size_t ld_domain = ((n_dofs_domain - 1) / align_elem + 1) * align_elem;

        // temporary allocations using the memory pool
        {
            cbmba<T> ator_dt(*cbmba_res_device, align_B);
            cbmba_res_device->do_transaction([&](){
                if(need_d_Us_dn) d_Us_dn[d] = std::make_unique<MatrixDense<T,I,cbmba>>(n_dofs_domain,    n_dofs_domain, ld_domain, true, ator_dt);
                if(need_d_Ls_dn) d_Ls_dn[d] = std::make_unique<MatrixDense<T,I,cbmba>>(n_dofs_domain,    n_dofs_domain, ld_domain, true, ator_dt);
                if(need_d_Xs)       d_Xs[d] = std::make_unique<MatrixDense<T,I,cbmba>>(n_dofs_interface, n_dofs_domain, ld_domain, true, ator_dt);
                if(need_d_Ys)       d_Ys[d] = std::make_unique<MatrixDense<T,I,cbmba>>(n_dofs_interface, n_dofs_domain, ld_domain, true, ator_dt);
                if(buffersizes[d] > 0) buffers[d] = cbmba_res_device->allocate(buffersizes[d], align_B);
            });
        }
        CHECK(rocblas_set_workspace(handles_blas[d], buffers[d], buffersizes[d]));

        // set the pointers inside the descriptors of the rest of the matrices
        {
            if(need_d_Us_dn) CHECK(rocsparse_dnmat_set_values(gpu_Us_dn[d], d_Us_dn[d]->vals));
            if(need_d_Ls_dn) CHECK(rocsparse_dnmat_set_values(gpu_Ls_dn[d], d_Ls_dn[d]->vals));
            if(need_d_Xs) CHECK(rocsparse_dnmat_set_values(gpu_Xs[d],    d_Xs[d]->vals));
            if(need_d_Xs) CHECK(rocsparse_dnmat_set_values(gpu_Xs_tc[d], d_Xs[d]->vals));
            if(need_d_Ys) CHECK(rocsparse_dnmat_set_values(gpu_Ys_tc[d], d_Ys[d]->vals));
            for(I j = 0; j < n_dofs_interface; j++) if(need_d_Xs) CHECK(rocsparse_dnvec_set_values(gpu_Xs_tc_cols[d][j], d_Xs[d]->vals + j * d_Xs[d]->ld));
            for(I j = 0; j < n_dofs_interface; j++) if(need_d_Ys) CHECK(rocsparse_dnvec_set_values(gpu_Ys_tc_cols[d][j], d_Ys[d]->vals + j * d_Ys[d]->ld));
        }

        // numeric factorization
        {
            cholmod_sparse cm_Kreg = my_cholmod_sparse_view(Kregs[d], 1);
            my_cholmod_factorize<I>(&cm_Kreg, cm_factors[d], cm_commons[d]);
        }

        // prepare matrices on host
        {
            my_cholmod_extract_factor_matrix(cm_factors[d], h_Us_sp[d], cm_commons[d]);
            if(need_h_Ls_sp) matrix_transpose(h_Ls_sp[d], h_Us_sp[d]);
        }

        // copy the new factors to device
        {
            if(need_d_Us_sp) copy_matrix_submit(d_Us_sp[d], h_Us_sp[d], hipMemcpyHostToDevice, streams[d]);
            if(need_d_Ls_sp) copy_matrix_submit(d_Ls_sp[d], h_Ls_sp[d], hipMemcpyHostToDevice, streams[d]);
            if(waitcontpar == 'W') CHECK(hipStreamSynchronize(streams[d]));
        }

        // prepare matrices on device
        {
            if(need_d_Us_dn) CHECK(rocsparse_sparse_to_dense(handles_sparse[d], gpu_Us_sp1[d],    gpu_Us_dn[d], rocsparse_sparse_to_dense_alg_default, &buffersizes[d], buffers[d]));
            if(need_d_Ls_dn) CHECK(rocsparse_sparse_to_dense(handles_sparse[d], gpu_Ls_sp[d],     gpu_Ls_dn[d], rocsparse_sparse_to_dense_alg_default, &buffersizes[d], buffers[d]));
            if(need_d_Xs)    CHECK(rocsparse_sparse_to_dense(handles_sparse[d], gpu_Bperms_sp[d], gpu_Xs[d],    rocsparse_sparse_to_dense_alg_default, &buffersizes[d], buffers[d]));
            if(waitcontpar == 'W') CHECK(hipStreamSynchronize(streams[d]));
        }

        // perform the actual assembly
        {
            T one = 1;
            T zero = 0;

            if(trsvtrsm == 'V' && spdnfactor == 'S' && transcallman == 'C') CHECK(rocsparse_spsv(handles_sparse[d], rocsparse_operation_transpose, &one, gpu_Us_sp1[d], gpu_Xs_tc_cols[d][0], gpu_Ys_tc_cols[d][0], rocsparse_data_type<T>(), rocsparse_spsv_alg_default, rocsparse_spsv_stage_preprocess, &buffersizes[d], buffers[d]));
            if(trsvtrsm == 'V' && spdnfactor == 'S' && transcallman == 'M') CHECK(rocsparse_spsv(handles_sparse[d], rocsparse_operation_none,      &one, gpu_Ls_sp[d],  gpu_Xs_tc_cols[d][0], gpu_Ys_tc_cols[d][0], rocsparse_data_type<T>(), rocsparse_spsv_alg_default, rocsparse_spsv_stage_preprocess, &buffersizes[d], buffers[d]));
            if(trsvtrsm == 'V' && spdnfactor == 'D') ;
            if(trsvtrsm == 'M' && spdnfactor == 'S' && transcallman == 'C') CHECK(rocsparse_spsm(handles_sparse[d], rocsparse_operation_transpose, rocsparse_operation_none, &one, gpu_Us_sp1[d], gpu_Xs_tc[d], gpu_Ys_tc[d], rocsparse_data_type<T>(), rocsparse_spsm_alg_default, rocsparse_spsm_stage_preprocess, &buffersizes[d], buffers[d]));
            if(trsvtrsm == 'M' && spdnfactor == 'S' && transcallman == 'M') CHECK(rocsparse_spsm(handles_sparse[d], rocsparse_operation_none,      rocsparse_operation_none, &one, gpu_Ls_sp[d],  gpu_Xs_tc[d], gpu_Ys_tc[d], rocsparse_data_type<T>(), rocsparse_spsm_alg_default, rocsparse_spsm_stage_preprocess, &buffersizes[d], buffers[d]));
            if(trsvtrsm == 'M' && spdnfactor == 'D') ;
            if(waitcontpar == 'W') CHECK(hipStreamSynchronize(streams[d]));

            if(trsvtrsm == 'V' && spdnfactor == 'S' && transcallman == 'C') for(I j = 0; j < n_dofs_interface; j++) CHECK(rocsparse_spsv(handles_sparse[d], rocsparse_operation_transpose, &one, gpu_Us_sp1[d], gpu_Xs_tc_cols[d][j], gpu_Ys_tc_cols[d][j], rocsparse_data_type<T>(), rocsparse_spsv_alg_default, rocsparse_spsv_stage_compute, &buffersizes[d], buffers[d]));
            if(trsvtrsm == 'V' && spdnfactor == 'S' && transcallman == 'M') for(I j = 0; j < n_dofs_interface; j++) CHECK(rocsparse_spsv(handles_sparse[d], rocsparse_operation_none,      &one, gpu_Ls_sp[d],  gpu_Xs_tc_cols[d][j], gpu_Ys_tc_cols[d][j], rocsparse_data_type<T>(), rocsparse_spsv_alg_default, rocsparse_spsv_stage_compute, &buffersizes[d], buffers[d]));
            if(trsvtrsm == 'V' && spdnfactor == 'D') for(I j = 0; j < n_dofs_interface; j++) CHECK(rocblas_dtrsv(handles_blas[d], rocblas_fill_lower, rocblas_operation_none, rocblas_diagonal_non_unit, n_dofs_domain, d_Us_dn[d]->vals, d_Us_dn[d]->ld, d_Xs[d]->vals + d_Xs[d]->ld * j, 1));
            if(trsvtrsm == 'M' && spdnfactor == 'S' && transcallman == 'C') CHECK(rocsparse_spsm(handles_sparse[d], rocsparse_operation_transpose, rocsparse_operation_none, &one, gpu_Us_sp1[d], gpu_Xs_tc[d], gpu_Ys_tc[d], rocsparse_data_type<T>(), rocsparse_spsm_alg_default, rocsparse_spsm_stage_compute, &buffersizes[d], buffers[d]));
            if(trsvtrsm == 'M' && spdnfactor == 'S' && transcallman == 'M') CHECK(rocsparse_spsm(handles_sparse[d], rocsparse_operation_none,      rocsparse_operation_none, &one, gpu_Ls_sp[d],  gpu_Xs_tc[d], gpu_Ys_tc[d], rocsparse_data_type<T>(), rocsparse_spsm_alg_default, rocsparse_spsm_stage_compute, &buffersizes[d], buffers[d]));
            if(trsvtrsm == 'M' && spdnfactor == 'D') CHECK(rocblas_dtrsm(handles_blas[d], rocblas_side_left, rocblas_fill_lower, rocblas_operation_none, rocblas_diagonal_non_unit, n_dofs_domain, n_dofs_interface, &one, d_Us_dn[d]->vals, d_Us_dn[d]->ld, d_Xs[d]->vals, d_Xs[d]->ld));
            if(waitcontpar == 'W') CHECK(hipStreamSynchronize(streams[d]));

            if(gemmsyrk == 'G')
            {
                if(trsvtrsm == 'V' && spdnfactor == 'S') CHECK(rocsparse_spsv(handles_sparse[d], rocsparse_operation_none, &one, gpu_Us_sp2[d], gpu_Ys_tc_cols[d][0], gpu_Xs_tc_cols[d][0], rocsparse_data_type<T>(), rocsparse_spsv_alg_default, rocsparse_spsv_stage_preprocess, &buffersizes[d], buffers[d]));
                if(trsvtrsm == 'V' && spdnfactor == 'D' && transcallman == 'C') ;
                if(trsvtrsm == 'V' && spdnfactor == 'D' && transcallman == 'M') ;
                if(trsvtrsm == 'M' && spdnfactor == 'S') CHECK(rocsparse_spsm(handles_sparse[d], rocsparse_operation_none, rocsparse_operation_none, &one, gpu_Us_sp2[d], gpu_Ys_tc[d], gpu_Xs_tc[d], rocsparse_data_type<T>(), rocsparse_spsm_alg_default, rocsparse_spsm_stage_preprocess, &buffersizes[d], buffers[d]));
                if(trsvtrsm == 'M' && spdnfactor == 'D' && transcallman == 'C') ;
                if(trsvtrsm == 'M' && spdnfactor == 'D' && transcallman == 'M') ;
                if(waitcontpar == 'W') CHECK(hipStreamSynchronize(streams[d]));

                if(trsvtrsm == 'V' && spdnfactor == 'S') for(I j = 0; j < n_dofs_interface; j++) CHECK(rocsparse_spsv(handles_sparse[d], rocsparse_operation_none, &one, gpu_Us_sp2[d], gpu_Ys_tc_cols[d][j], gpu_Xs_tc_cols[d][j], rocsparse_data_type<T>(), rocsparse_spsv_alg_default, rocsparse_spsv_stage_compute, &buffersizes[d], buffers[d]));
                if(trsvtrsm == 'V' && spdnfactor == 'D' && transcallman == 'C') for(I j = 0; j < n_dofs_interface; j++) CHECK(rocblas_dtrsv(handles_blas[d], rocblas_fill_lower, rocblas_operation_transpose, rocblas_diagonal_non_unit, n_dofs_domain, d_Us_dn[d]->vals, d_Us_dn[d]->ld, d_Xs[d]->vals + d_Xs[d]->ld * j, 1));
                if(trsvtrsm == 'V' && spdnfactor == 'D' && transcallman == 'M') for(I j = 0; j < n_dofs_interface; j++) CHECK(rocblas_dtrsv(handles_blas[d], rocblas_fill_upper, rocblas_operation_none,      rocblas_diagonal_non_unit, n_dofs_domain, d_Ls_dn[d]->vals, d_Ls_dn[d]->ld, d_Xs[d]->vals + d_Xs[d]->ld * j, 1));
                if(trsvtrsm == 'M' && spdnfactor == 'S') CHECK(rocsparse_spsm(handles_sparse[d], rocsparse_operation_none, rocsparse_operation_none, &one, gpu_Us_sp2[d], gpu_Ys_tc[d], gpu_Xs_tc[d], rocsparse_data_type<T>(), rocsparse_spsm_alg_default, rocsparse_spsm_stage_compute, &buffersizes[d], buffers[d]));
                if(trsvtrsm == 'M' && spdnfactor == 'D' && transcallman == 'C') CHECK(rocblas_dtrsm(handles_blas[d], rocblas_side_left, rocblas_fill_lower, rocblas_operation_transpose, rocblas_diagonal_non_unit, n_dofs_domain, n_dofs_interface, &one, d_Us_dn[d]->vals, d_Us_dn[d]->ld, d_Xs[d]->vals, d_Xs[d]->ld));
                if(trsvtrsm == 'M' && spdnfactor == 'D' && transcallman == 'M') CHECK(rocblas_dtrsm(handles_blas[d], rocblas_side_left, rocblas_fill_upper, rocblas_operation_none,      rocblas_diagonal_non_unit, n_dofs_domain, n_dofs_interface, &one, d_Ls_dn[d]->vals, d_Ls_dn[d]->ld, d_Xs[d]->vals, d_Xs[d]->ld));
                if(waitcontpar == 'W') CHECK(hipStreamSynchronize(streams[d]));

                CHECK(rocsparse_spmm(handles_sparse[d], rocsparse_operation_none, rocsparse_operation_none, &one, gpu_Bperms_sp[d], gpu_Xs_tc[d], &zero, gpu_Fs[d], rocsparse_data_type<T>(), rocsparse_spmm_alg_default, rocsparse_spmm_stage_preprocess, &buffersizes[d], buffers[d]));
                CHECK(rocsparse_spmm(handles_sparse[d], rocsparse_operation_none, rocsparse_operation_none, &one, gpu_Bperms_sp[d], gpu_Xs_tc[d], &zero, gpu_Fs[d], rocsparse_data_type<T>(), rocsparse_spmm_alg_default, rocsparse_spmm_stage_compute,    &buffersizes[d], buffers[d]));
                if(waitcontpar == 'W') CHECK(hipStreamSynchronize(streams[d]));
            }
            if(gemmsyrk == 'S')
            {
                if(spdnfactor == 'S') CHECK(rocblas_dsyrk(handles_blas[d], rocblas_fill_lower, rocblas_operation_transpose, n_dofs_interface, n_dofs_domain, &one, d_Ys[d]->vals, d_Ys[d]->ld, &zero, d_Fs[d].vals, d_Fs[d].ld));
                if(spdnfactor == 'D') CHECK(rocblas_dsyrk(handles_blas[d], rocblas_fill_lower, rocblas_operation_transpose, n_dofs_interface, n_dofs_domain, &one, d_Xs[d]->vals, d_Xs[d]->ld, &zero, d_Fs[d].vals, d_Fs[d].ld));
                if(waitcontpar == 'W') CHECK(hipStreamSynchronize(streams[d]));
            }
        }

        // free the temporary memory from the pool
        my_hip_submit_host_function(streams[d], [&,d](){
            if(need_d_Us_dn) { d_Us_dn[d]->deallocate(); d_Us_dn[d] = nullptr; };
            if(need_d_Ls_dn) { d_Ls_dn[d]->deallocate(); d_Ls_dn[d] = nullptr; };
            if(need_d_Xs)    {    d_Xs[d]->deallocate();    d_Xs[d] = nullptr; };
            if(need_d_Ys)    {    d_Ys[d]->deallocate();    d_Ys[d] = nullptr; };
            cbmba_res_device->deallocate(buffers[d]); buffers[d] = nullptr;
        });
        if(waitcontpar == 'W') CHECK(hipStreamSynchronize(streams[d]));
    }

    // clean up the mess from buggy openmp in clang
    run_dummy_parallel_region();

    CHECK(hipDeviceSynchronize());

    stage = 2;
}


template<typename T, typename I>
void my_dual_operator_cluster<T,I>::apply(MatrixDense<T,I> & y_cluster, const MatrixDense<T,I> & x_cluster, const std::vector<std::vector<I>> & domain_to_cluster_maps)
{
    // y = F * x
    // temporary solution, will be better

    char waitcontpar  = ms_get_waitcontpar (magicstring);

    if(stage != 2) throw std::runtime_error("invalid stage when calling apply");

    #pragma omp parallel for schedule(static,1) if(waitcontpar == 'P')
    for(size_t d = 0; d < n_domains; d++)
    {
        T one = 1;
        T zero = 0;
        I n_dofs_interface = h_Bperms_sp[d].nrows;

        for(I i = 0; i < n_dofs_interface; i++) h_apply_xs[d].vals[i] = x_cluster.vals[domain_to_cluster_maps[d][i]];

        copy_matrix_submit(d_apply_xs[d], h_apply_xs[d], streams[d]);
        if(waitcontpar == 'W') CHECK(hipStreamSynchronize(streams[d]));

        CHECK(rocblas_dsymv(handles_blas[d], rocblas_fill_lower, d_Fs[d].nrows, &one, d_Fs[d].vals, d_Fs[d].ld, d_apply_xs[d].vals, 1, &zero, d_apply_ys[d].vals, 1));
        if(waitcontpar == 'W') CHECK(hipStreamSynchronize(streams[d]));

        copy_matrix_submit(h_apply_ys[d], d_apply_ys[d], streams[d]);
        if(waitcontpar == 'W') CHECK(hipStreamSynchronize(streams[d]));
    }

    for(I i = 0; i < y_cluster.nrows; i++) y_cluster.vals[i] = T{0};

    CHECK(hipDeviceSynchronize());

    #pragma omp parallel for schedule(static,1) if(waitcontpar == 'P')
    for(size_t d = 0; d < n_domains; d++)
    {
        I n_dofs_interface = h_Bperms_sp[d].nrows;
        for(I i = 0; i < n_dofs_interface; i++)
        {
            #pragma omp atomic
            y_cluster.vals[domain_to_cluster_maps[d][i]] += h_apply_ys[d].vals[i];
        }
    }
}


template<typename T, typename I>
void my_dual_operator_cluster<T,I>::destroy()
{
    if(stage != 1 && stage != 2) throw std::runtime_error("invalid stage when calling destroy");

    for(size_t d = 0; d < n_domains; d++)
    {
        I n_dofs_interface = h_Bperms_sp[d].nrows;

        CHECK(rocsparse_destroy_spmat_descr(gpu_Us_sp1[d]));
        CHECK(rocsparse_destroy_spmat_descr(gpu_Us_sp2[d]));
        CHECK(rocsparse_destroy_spmat_descr(gpu_Ls_sp[d]));
        CHECK(rocsparse_destroy_spmat_descr(gpu_Bperms_sp[d]));
        CHECK(rocsparse_destroy_dnmat_descr(gpu_Us_dn[d]));
        CHECK(rocsparse_destroy_dnmat_descr(gpu_Ls_dn[d]));
        CHECK(rocsparse_destroy_dnmat_descr(gpu_Xs_tc[d]));
        CHECK(rocsparse_destroy_dnmat_descr(gpu_Ys_tc[d]));
        for(I j = 0; j < n_dofs_interface; j++) CHECK(rocsparse_destroy_dnvec_descr(gpu_Xs_tc_cols[d][j]));
        for(I j = 0; j < n_dofs_interface; j++) CHECK(rocsparse_destroy_dnvec_descr(gpu_Ys_tc_cols[d][j]));
        CHECK(rocsparse_destroy_dnmat_descr(gpu_Fs[d]));

        my_cholmod_free_factor<I>(&cm_factors[d], cm_commons[d]);
        my_cholmod_finish_and_delete(&cm_commons[d]);

        CHECK(rocsparse_destroy_handle(handles_sparse[d]));
        CHECK(rocblas_destroy_handle(handles_blas[d]));
        CHECK(hipStreamDestroy(streams[d]));
    }
    CHECK(hipFree(mem_pool_device));

    mem_pool_device = nullptr;
    cbmba_res_device = nullptr;
    streams.clear();
    handles_sparse.clear();
    handles_blas.clear();
    d_Fs.clear();
    h_Us_sp.clear();
    h_Ls_sp.clear();
    h_Bperms_sp.clear();
    d_Us_sp.clear();
    d_Ls_sp.clear();
    d_Bperms_sp.clear();
    d_Us_dn.clear();
    d_Ls_dn.clear();
    d_Xs.clear();
    d_Ys.clear();
    gpu_Fs.clear();
    gpu_Us_sp1.clear();
    gpu_Us_sp2.clear();
    gpu_Ls_sp.clear();
    gpu_Bperms_sp.clear();
    gpu_Us_dn.clear();
    gpu_Ls_dn.clear();
    gpu_Xs.clear();
    gpu_Xs_tc.clear();
    gpu_Ys_tc.clear();
    gpu_Xs_tc_cols.clear();
    gpu_Ys_tc_cols.clear();
    cm_commons.clear();
    cm_factors.clear();
    buffers.clear();
    buffersizes.clear();
    d_apply_xs.clear();
    d_apply_ys.clear();
    h_apply_xs.clear();
    h_apply_ys.clear();

    CHECK(hipDeviceSynchronize());

    stage = 0;
}

/*
 * void set(const std::vector<MatrixCSR<T,I>> & Kregs, const std::vector<MatrixCSR<T,I>> & Bs);
   void update(const std::vector<MatrixCSR<T,I>> & Kregs);
   void apply(MatrixDense<T,I> & y_cluster, const MatrixDense<T,I> & x_cluster, const std::vector<std::vector<I>> & domain_to_cluster_maps);
   void destroy();
 */

template <typename T>
struct Acc_FETI_Dual_Operator {
	my_dual_operator_cluster<T, int> dual;

	Acc_FETI_Dual_Operator(): dual("___PPMS_MMS_") {}
};

template <>
AccFETIDualOperator<double, Matrix_CSR>::AccFETIDualOperator()
: _acc(nullptr)
{
	_acc = new Acc_FETI_Dual_Operator<double>();
}


template <>
AccFETIDualOperator<double, Matrix_CSR>::~AccFETIDualOperator()
{
	_acc->dual.destroy();
}

template <>
void AccFETIDualOperator<double, Matrix_CSR>::set(const std::vector<Matrix_CSR<double> > &K, const std::vector<Matrix_CSR<double> > &B)
{
	std::vector<MatrixCSR<double,int> > _K, _B;
	_K.resize(K.size());
	_B.resize(B.size());
	for (size_t di = 0; di < K.size(); ++di) {
		printf("%lu %d %d %d\n", di, K[di].nrows, K[di].ncols, K[di].nnz);
		_K[di].resize(K[di].nrows, K[di].ncols, K[di].nnz, false);
		_K[di].rowptrs = K[di].rows;
		_K[di].colidxs = K[di].cols;
		_K[di].vals = K[di].vals;

		printf("%lu %d %d %d\n", di, B[di].nrows, B[di].ncols, B[di].nnz);
		_B[di].resize(B[di].nrows, B[di].ncols, B[di].nnz, false);
		_B[di].rowptrs = B[di].rows;
		_B[di].colidxs = B[di].cols;
		_B[di].vals = B[di].vals;
	}
	_acc->dual.set(_K, _B);
}

template <>
void AccFETIDualOperator<double, Matrix_CSR>::update(const std::vector<Matrix_CSR<double> > &K)
{
	std::vector<MatrixCSR<double,int> > _K;
	_K.resize(K.size());
	for (size_t di = 0; di < K.size(); ++di) {
		_K[di].resize(K[di].nrows, K[di].ncols, K[di].nnz, false);
		_K[di].rowptrs = K[di].rows;
		_K[di].colidxs = K[di].cols;
		_K[di].vals = K[di].vals;
	}

	_acc->dual.update(_K);
}

template <>
void AccFETIDualOperator<double, Matrix_CSR>::apply(const Vector_Dual<double> &x, Vector_Dual<double> &y, const std::vector<std::vector<int> > & D2C)
{
	MatrixDense<double, int> _x, _y;
	_x.resize(x.size, 1, -1, false);
	_x.vals = x.vals;
	_y.resize(y.size, 1, -1, false);
	_y.vals = y.vals;
	_acc->dual.apply(_x, _y, D2C);
}

template <typename T, template <typename> class Matrix>
AccFETIDualOperator<T, Matrix>::AccFETIDualOperator()
: _acc(nullptr)
{
}


template <typename T, template <typename> class Matrix>
AccFETIDualOperator<T, Matrix>::~AccFETIDualOperator()
{
}

template <typename T, template <typename> class Matrix>
void AccFETIDualOperator<T, Matrix>::set(const std::vector<Matrix<T> > &K, const std::vector<Matrix<T> > &B)
{
}

template <typename T, template <typename> class Matrix>
void AccFETIDualOperator<T, Matrix>::update(const std::vector<Matrix<T> > &K)
{
}

template <typename T, template <typename> class Matrix>
void AccFETIDualOperator<T, Matrix>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y, const std::vector<std::vector<int> > & D2C)
{
}

template struct AccFETIDualOperator<double, Matrix_CSR>;
template struct AccFETIDualOperator<std::complex<double>, Matrix_CSR>;

}


#endif
