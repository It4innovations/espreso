
#include "totalfeti.explicit.acc.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/mpiinfo.h"

#include "my_timer.h"

namespace espreso {

template <typename T, typename I>
TotalFETIExplicitAcc<T,I>::TotalFETIExplicitAcc(FETI<T> &feti)
: DualOperator<T>(feti), n_domains(0), n_queues(0), mem_pool_device(nullptr)
{
    if(stage != 0) eslog::error("init: invalid order of operations in dualop\n");

    device = gpu::mgm::get_device_by_mpi(info::mpi::rank, info::mpi::size);

    const char * magicstring = "PPCSRLS_TG";

    wcpset       = magicstring[0];
    wcpupdate    = magicstring[1];
    wcpapply     = magicstring[2];
    spdnfactor   = magicstring[3];
    trsvtrsm     = magicstring[4];
    factrs1      = magicstring[5];
    factrs2syrk  = magicstring[6];
    handlecount  = magicstring[8];
    applyalg     = magicstring[9];

    stage = 1;
}

template <typename T, typename I>
TotalFETIExplicitAcc<T,I>::~TotalFETIExplicitAcc()
{
    gpu::mgm::set_device(device);

    my_timer tm_total, tm_descriptors, tm_gpulibs, tm_gpumempool, tm_gpumem, tm_gpuhostmem, tm_wait;

    tm_total.start();

    tm_descriptors.start();
    for(size_t d = 0; d < n_domains; d++)
    {
        I n_dofs_interface = h_Bperms_sp[d].nrows;

        if(!descr_Us_sp1.empty())    gpu::spblas::descr_matrix_csr_destroy(descr_Us_sp1[d]);
        if(!descr_Us_sp2.empty())    gpu::spblas::descr_matrix_csr_destroy(descr_Us_sp2[d]);
        if(!descr_Ls_sp1.empty())    gpu::spblas::descr_matrix_csr_destroy(descr_Ls_sp1[d]);
        if(!descr_Ls_sp2.empty())    gpu::spblas::descr_matrix_csr_destroy(descr_Ls_sp2[d]);
        if(!descr_Bperms_sp.empty()) gpu::spblas::descr_matrix_csr_destroy(descr_Bperms_sp[d]);
        if(!descr_Us_dn.empty()) gpu::spblas::descr_matrix_dense_destroy(descr_Us_dn[d]);
        if(!descr_Ls_dn.empty()) gpu::spblas::descr_matrix_dense_destroy(descr_Ls_dn[d]);
        if(!descr_Xs_c.empty()) gpu::spblas::descr_matrix_dense_destroy(descr_Xs_c[d]);
        if(!descr_Xs_r.empty()) gpu::spblas::descr_matrix_dense_destroy(descr_Xs_r[d]);
        if(!descr_Ys_c.empty()) gpu::spblas::descr_matrix_dense_destroy(descr_Ys_c[d]);
        if(!descr_Ys_r.empty()) gpu::spblas::descr_matrix_dense_destroy(descr_Ys_r[d]);
        if(!descr_Xs_vecs.empty()) for(I j = 0; j < n_dofs_interface; j++) gpu::spblas::descr_vector_dense_destroy(descr_Xs_vecs[d][j]);
        if(!descr_Ys_vecs.empty()) for(I j = 0; j < n_dofs_interface; j++) gpu::spblas::descr_vector_dense_destroy(descr_Ys_vecs[d][j]);
        if(!descr_Fs_c.empty()) gpu::spblas::descr_matrix_dense_destroy(descr_Fs_c[d]);
        if(!descr_Fs_r.empty()) gpu::spblas::descr_matrix_dense_destroy(descr_Fs_r[d]);
        if(!descrs_sparse_trsv1.empty()) gpu::spblas::descr_sparse_trsv_destroy(descrs_sparse_trsv1[d]);
        if(!descrs_sparse_trsm1.empty()) gpu::spblas::descr_sparse_trsm_destroy(descrs_sparse_trsm1[d]);
        if(!descrs_sparse_trsv2.empty()) gpu::spblas::descr_sparse_trsv_destroy(descrs_sparse_trsv2[d]);
        if(!descrs_sparse_trsm2.empty()) gpu::spblas::descr_sparse_trsm_destroy(descrs_sparse_trsm2[d]);
    }
    tm_descriptors.stop();

    tm_gpulibs.start();
    gpu::mgm::queue_destroy(main_q);
    for(gpu::mgm::queue & q : queues) gpu::mgm::queue_destroy(q);
    for(gpu::dnblas::handle & h : handles_dense) gpu::dnblas::handle_destroy(h);
    for(gpu::spblas::handle & h : handles_sparse) gpu::spblas::handle_destroy(h);
    tm_gpulibs.stop();

    tm_gpumempool.start();
    gpu::mgm::memfree_device(mem_pool_device);
    tm_gpumempool.stop();

    tm_gpumem.start();
    d_Us_sp.clear();
    d_Ls_sp.clear();
    d_Bperms_sp.clear();
    d_Us_dn.clear();
    d_Ls_dn.clear();
    d_Xs_c.clear();
    d_Xs_r.clear();
    d_Ys_c.clear();
    d_Ys_r.clear();
    d_Fs.clear();
    for(void * ptr : buffers_sptrs1) gpu::mgm::memfree_device(ptr);
    for(void * ptr : buffers_sptrs2) gpu::mgm::memfree_device(ptr);
    for(void * ptr : buffers_spmm) gpu::mgm::memfree_device(ptr);
    d_applyg_D2Cs.clear();
    d_apply_xs.clear();
    d_apply_ys.clear();
    d_applyg_x_cluster.clear();
    d_applyg_y_cluster.clear();
    d_applyg_xs_pointers.clear();
    d_applyg_ys_pointers.clear();
    d_applyg_n_dofs_interfaces.clear();
    d_applyg_D2Cs_pointers.clear();
    tm_gpumem.stop();

    tm_gpuhostmem.start();
    h_Us_sp.clear();
    h_Ls_sp.clear();
    h_Bperms_sp.clear();
    h_applyc_xs.clear();
    h_applyc_ys.clear();
    tm_gpuhostmem.stop();

    tm_wait.start();
    gpu::mgm::device_wait();
    tm_wait.stop();

    tm_total.stop();

    stage = 0;

    print_timer("Destroy total", tm_total);
    print_timer("Destroy   descriptors", tm_descriptors);
    print_timer("Destroy   gpulibs", tm_gpulibs);
    print_timer("Destroy   gpumempool", tm_gpumempool);
    print_timer("Destroy   gpumem", tm_gpumem);
    print_timer("Destroy   gpuhostmem", tm_gpuhostmem);
    print_timer("Destroy   wait", tm_wait);
}

template <typename T, typename I>
void TotalFETIExplicitAcc<T,I>::info()
{
    // DualOperatorInfo sum, min, max;
    size_t minF = INT32_MAX, maxF = 0, sumF = 0;
    for (size_t d = 0; d < d_Fs.size(); ++d) {
        minF = std::min(minF, d_Fs[d].nrows * d_Fs[d].ncols * sizeof(T));
        maxF = std::max(maxF, d_Fs[d].nrows * d_Fs[d].ncols * sizeof(T));
        sumF += d_Fs[d].nrows * d_Fs[d].ncols * sizeof(T);
    }

//    TotalFETIImplicit<T>::reduceInfo(sum, min, max);
    Communication::allReduce(&minF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_MIN);
    Communication::allReduce(&maxF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_MAX);
    Communication::allReduce(&sumF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_SUM);

    eslog::info(" = EXPLICIT TOTAL FETI OPERATOR ON GPU                                                       = \n");
//    TotalFETIImplicit<T>::printInfo(sum, min, max);
    eslog::info(" =   F MEMORY [MB]                                            %8.2f <%8.2f - %8.2f> = \n", (double)sumF / d_Fs.size() / 1024. / 1024., minF / 1024. / 1024., maxF / 1024. / 1024.);
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}





template <typename T, typename I>
void TotalFETIExplicitAcc<T,I>::set(const step::Step &step)
{
    if(stage != 1) eslog::error("set: invalid order of operations in dualop\n");
    if(DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::NONSYMMETRIC_BOTH && (factrs1 != 'L' || factrs2syrk != 'U')) eslog::error("Wrong options for non-symmetric K\n");

    my_timer tm_total, tm_gpuinit, tm_gpuset, tm_vecresize, tm_gpucreate, tm_mainloop_outer, tm_mainloop_inner, tm_Kreg_combine, tm_solver_commit, tm_fact_symbolic, tm_descriptors, tm_buffersize, tm_alloc, tm_alloc_host, tm_alloc_device, tm_setpointers, tm_Bperm, tm_get_factors, tm_extract, tm_transpose, tm_copyin, tm_kernels_preprocess, tm_trs1, tm_trs2, tm_gemm, tm_applystuff, tm_poolalloc, tm_wait;

    tm_total.start();
    n_domains = feti.K.size();

    if(handlecount == 'D') n_queues = n_domains;
    if(handlecount == 'T') n_queues = omp_get_max_threads();

    const bool need_L = (factrs1 == 'L' || factrs2syrk == 'L');
    const bool need_U = (factrs1 == 'U' || factrs2syrk == 'U');
    const bool need_d_Ls_dn = (need_L && spdnfactor == 'D');
    const bool need_d_Us_dn = (need_U && spdnfactor == 'D');
    const bool need_d_Ls_sp = need_L;
    const bool need_d_Us_sp = need_U;
    const bool need_h_Ls_sp = (need_d_Ls_sp || DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_LOWER || DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::NONSYMMETRIC_BOTH);
    const bool need_h_Us_sp = (need_d_Us_sp || DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_UPPER || DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::NONSYMMETRIC_BOTH);
    const bool need_d_Bperms_sp = true;
    const bool need_h_Bperms_sp = true;
    const bool need_d_Xs_c = (trsvtrsm == 'C' || trsvtrsm == 'V');
    const bool need_d_Xs_r = (trsvtrsm == 'R');
    const bool need_d_Ys_c = (need_d_Xs_c && spdnfactor == 'S');
    const bool need_d_Ys_r = (need_d_Xs_r && spdnfactor == 'S');
    const bool need_d_Xs_vecs = (trsvtrsm == 'V');
    const bool need_d_Ys_vecs = (need_d_Xs_vecs && spdnfactor == 'S');

    tm_gpuinit.start();
    gpu::mgm::init_gpu(device);
    tm_gpuinit.stop();

    tm_gpuset.start();
    gpu::mgm::set_device(device);
    tm_gpuset.stop();

    // approximate check if it is possible for the matrices to fit into memory
    {
        size_t mem_needed = 0;
        for(size_t d = 0; d < n_domains; d++)
        {
            I n_dofs_interface = feti.B1[d].nrows;
            I ld_interface = ((n_dofs_interface - 1) / align_elem + 1) * align_elem;
            Matrix_Dense<T,I>::memoryRequirement(n_dofs_interface, n_dofs_interface, ld_interface);
        }
        size_t mem_capacity = gpu::mgm::get_device_memory_capacity();
        if(mem_needed > mem_capacity)
        {
            eslog::error("Not enough device memory for this input. Need %zu GiB, but have only %zu GiB\n", mem_needed >> 30, mem_capacity >> 30);
        }
    }

    tm_vecresize.start();
    queues.resize(n_queues);
    handles_dense.resize(n_queues);
    handles_sparse.resize(n_queues);
    solvers_Kreg.resize(n_domains);
    Kregs.resize(n_domains);
    d_Fs.resize(n_domains);
    if(need_h_Us_sp) h_Us_sp.resize(n_domains);
    if(need_h_Ls_sp) h_Ls_sp.resize(n_domains);
    if(need_h_Bperms_sp) h_Bperms_sp.resize(n_domains);
    if(need_d_Us_sp) d_Us_sp.resize(n_domains);
    if(need_d_Ls_sp) d_Ls_sp.resize(n_domains);
    if(need_d_Bperms_sp) d_Bperms_sp.resize(n_domains);
    if(need_d_Us_dn) d_Us_dn.resize(n_domains);
    if(need_d_Ls_dn) d_Ls_dn.resize(n_domains);
    if(need_d_Xs_c) d_Xs_c.resize(n_domains);
    if(need_d_Xs_r) d_Xs_r.resize(n_domains);
    if(need_d_Ys_c) d_Ys_c.resize(n_domains);
    if(need_d_Ys_r) d_Ys_r.resize(n_domains);
    if(true) descr_Fs_r.resize(n_domains);
    if(true) descr_Fs_c.resize(n_domains);
    if(need_d_Us_sp) descr_Us_sp1.resize(n_domains);
    if(need_d_Us_sp) descr_Us_sp2.resize(n_domains);
    if(need_d_Ls_sp) descr_Ls_sp1.resize(n_domains);
    if(need_d_Ls_sp) descr_Ls_sp2.resize(n_domains);
    if(need_d_Bperms_sp) descr_Bperms_sp.resize(n_domains);
    if(need_d_Us_dn) descr_Us_dn.resize(n_domains);
    if(need_d_Ls_dn) descr_Ls_dn.resize(n_domains);
    if(need_d_Xs_r) descr_Xs_r.resize(n_domains);
    if(need_d_Xs_c) descr_Xs_c.resize(n_domains);
    if(need_d_Ys_r) descr_Ys_r.resize(n_domains);
    if(need_d_Ys_c) descr_Ys_c.resize(n_domains);
    if(need_d_Xs_vecs) descr_Xs_vecs.resize(n_domains);
    if(need_d_Ys_vecs) descr_Ys_vecs.resize(n_domains);
    if(spdnfactor == 'S' && trsvtrsm == 'V') descrs_sparse_trsv1.resize(n_domains);
    if(spdnfactor == 'S' && trsvtrsm != 'V') descrs_sparse_trsm1.resize(n_domains);
    if(spdnfactor == 'S' && trsvtrsm == 'V' && factrs2syrk != 'S') descrs_sparse_trsv2.resize(n_domains);
    if(spdnfactor == 'S' && trsvtrsm != 'V' && factrs2syrk != 'S') descrs_sparse_trsm2.resize(n_domains);
    if(spdnfactor == 'S') buffersizes_sptrs1.resize(n_domains, 0);
    if(spdnfactor == 'S' && factrs2syrk != 'S') buffersizes_sptrs2.resize(n_domains, 0);
    if(factrs2syrk != 'S') buffersizes_spmm.resize(n_domains, 0);
    buffersizes_other.resize(n_domains, 0);
    if(spdnfactor == 'S') buffers_sptrs1.resize(n_domains, nullptr);
    if(spdnfactor == 'S' && factrs2syrk != 'S') buffers_sptrs2.resize(n_domains, nullptr);
    if(factrs2syrk != 'S') buffers_spmm.resize(n_domains, nullptr);
    transmaps_L2U.resize(n_domains);
    transmaps_U2L.resize(n_domains);
    d_applyg_D2Cs.resize(n_domains);
    d_apply_xs.resize(n_domains);
    d_apply_ys.resize(n_domains);
    h_applyc_xs.resize(n_domains);
    h_applyc_ys.resize(n_domains);
    tm_vecresize.stop();

    // create gpu stream and libraries
    tm_gpucreate.start();
    gpu::mgm::queue_create(main_q);
    for(gpu::mgm::queue & q : queues) gpu::mgm::queue_create(q);
    for(size_t i = 0; i < n_queues; i++) gpu::dnblas::handle_create(handles_dense[i], queues[i]);
    for(size_t i = 0; i < n_queues; i++) gpu::spblas::handle_create(handles_sparse[i], queues[i]);
    tm_gpucreate.stop();

    tm_mainloop_outer.start();
    #pragma omp parallel for schedule(static,1) if(wcpset == 'P')
    for(size_t d = 0; d < n_domains; d++)
    {
        tm_mainloop_inner.start();

        I n_dofs_domain = feti.B1[d].ncols;
        I n_dofs_interface = feti.B1[d].nrows;
        I ld_domain = ((n_dofs_domain - 1) / align_elem + 1) * align_elem;
        I ld_interface = ((n_dofs_interface - 1) / align_elem + 1) * align_elem;
        I n_nz_factor;

        gpu::mgm::queue & q = queues[d % n_queues];
        gpu::dnblas::handle & hd = handles_dense[d % n_queues];
        gpu::spblas::handle & hs = handles_sparse[d % n_queues];

        // Kreg = K + RegMat symbolic pattern
        tm_Kreg_combine.start();
        {
            math::combine(Kregs[d], feti.K[d], feti.RegMat[d]);
            if constexpr(utils::is_real<T>())    Kregs[d].type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
            if constexpr(utils::is_complex<T>()) Kregs[d].type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE;
            Kregs[d].shape = feti.K[d].shape;
        }
        tm_Kreg_combine.stop();

        // commit Kreg to solver (just symbolic pattern present now)
        tm_solver_commit.start();
        {
            solvers_Kreg[d].commit(Kregs[d]);
        }
        tm_solver_commit.stop();

        // symbolic factorization
        tm_fact_symbolic.start();
        {
            solvers_Kreg[d].symbolicFactorization();
            n_nz_factor = solvers_Kreg[d].getFactorNnz();
        }
        tm_fact_symbolic.stop();

        // create descriptors
        tm_descriptors.start();
        {
            if(need_d_Us_sp)     gpu::spblas::descr_matrix_csr_create<T,I>(   descr_Us_sp1[d], n_dofs_domain,    n_dofs_domain, n_nz_factor, 'U');
            if(need_d_Us_sp)     gpu::spblas::descr_matrix_csr_create<T,I>(   descr_Us_sp2[d], n_dofs_domain,    n_dofs_domain, n_nz_factor, 'U');
            if(need_d_Ls_sp)     gpu::spblas::descr_matrix_csr_create<T,I>(   descr_Ls_sp1[d], n_dofs_domain,    n_dofs_domain, n_nz_factor, 'L');
            if(need_d_Ls_sp)     gpu::spblas::descr_matrix_csr_create<T,I>(   descr_Ls_sp2[d], n_dofs_domain,    n_dofs_domain, n_nz_factor, 'L');
            if(need_d_Bperms_sp) gpu::spblas::descr_matrix_csr_create<T,I>(descr_Bperms_sp[d], n_dofs_interface, n_dofs_domain, feti.B1[d].nnz, 'N');
            if(need_d_Us_dn) gpu::spblas::descr_matrix_dense_create<T,I>(descr_Us_dn[d], n_dofs_domain,    n_dofs_domain,    ld_domain,    'R');
            if(need_d_Ls_dn) gpu::spblas::descr_matrix_dense_create<T,I>(descr_Ls_dn[d], n_dofs_domain,    n_dofs_domain,    ld_domain,    'R');
            if(need_d_Xs_c)  gpu::spblas::descr_matrix_dense_create<T,I>(descr_Xs_c[d], n_dofs_domain,    n_dofs_interface, ld_domain,    'C');
            if(need_d_Xs_r)  gpu::spblas::descr_matrix_dense_create<T,I>(descr_Xs_r[d], n_dofs_domain,    n_dofs_interface, ld_interface, 'R');
            if(need_d_Ys_c)  gpu::spblas::descr_matrix_dense_create<T,I>(descr_Ys_c[d], n_dofs_domain,    n_dofs_interface, ld_domain,    'C');
            if(need_d_Ys_r)  gpu::spblas::descr_matrix_dense_create<T,I>(descr_Ys_r[d], n_dofs_domain,    n_dofs_interface, ld_interface, 'R');
            if(need_d_Xs_vecs) descr_Xs_vecs[d].resize(n_dofs_interface);
            if(need_d_Xs_vecs) for(I j = 0; j < n_dofs_interface; j++) gpu::spblas::descr_vector_dense_create<T,I>(descr_Xs_vecs[d][j], n_dofs_domain);
            if(need_d_Ys_vecs) descr_Ys_vecs[d].resize(n_dofs_interface);
            if(need_d_Ys_vecs) for(I j = 0; j < n_dofs_interface; j++) gpu::spblas::descr_vector_dense_create<T,I>(descr_Ys_vecs[d][j], n_dofs_domain);
            if(true) gpu::spblas::descr_matrix_dense_create<T,I>(descr_Fs_r[d], n_dofs_interface, n_dofs_interface, ld_interface, 'R');
            if(true) gpu::spblas::descr_matrix_dense_create<T,I>(descr_Fs_c[d], n_dofs_interface, n_dofs_interface, ld_interface, 'C');
            if(spdnfactor == 'S' && trsvtrsm == 'V') gpu::spblas::descr_sparse_trsv_create(descrs_sparse_trsv1[d]);
            if(spdnfactor == 'S' && trsvtrsm != 'V') gpu::spblas::descr_sparse_trsm_create(descrs_sparse_trsm1[d]);
            if(spdnfactor == 'S' && trsvtrsm == 'V' && factrs2syrk != 'S') gpu::spblas::descr_sparse_trsv_create(descrs_sparse_trsv2[d]);
            if(spdnfactor == 'S' && trsvtrsm != 'V' && factrs2syrk != 'S') gpu::spblas::descr_sparse_trsm_create(descrs_sparse_trsm2[d]);
        }
        tm_descriptors.stop();

        // buffersize
        tm_buffersize.start();
        {
            std::vector<size_t> buffer_requirements;

            if(need_d_Us_dn) gpu::spblas::sparse_to_dense<T,I>(hs, 'N', descr_Us_sp1[d],    descr_Us_dn[d], buffer_requirements.emplace_back(), nullptr, 'B');
            if(need_d_Ls_dn) gpu::spblas::sparse_to_dense<T,I>(hs, 'N', descr_Ls_sp1[d],    descr_Ls_dn[d], buffer_requirements.emplace_back(), nullptr, 'B');
            if(need_d_Xs_c)  gpu::spblas::sparse_to_dense<T,I>(hs, 'T', descr_Bperms_sp[d], descr_Xs_c[d], buffer_requirements.emplace_back(), nullptr, 'B');
            if(need_d_Xs_r)  gpu::spblas::sparse_to_dense<T,I>(hs, 'T', descr_Bperms_sp[d], descr_Xs_r[d], buffer_requirements.emplace_back(), nullptr, 'B');

            if(trsvtrsm == 'V' && spdnfactor == 'S' && factrs1 == 'L') gpu::spblas::trsv<T,I>(hs, 'N', descr_Ls_sp1[d], descr_Xs_vecs[d][0], descr_Ys_vecs[d][0], descrs_sparse_trsv1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'B');
            if(trsvtrsm == 'V' && spdnfactor == 'S' && factrs1 == 'U') gpu::spblas::trsv<T,I>(hs, 'H', descr_Us_sp1[d], descr_Xs_vecs[d][0], descr_Ys_vecs[d][0], descrs_sparse_trsv1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'B');
            if(trsvtrsm == 'C' && spdnfactor == 'S' && factrs1 == 'L') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Ls_sp1[d], descr_Xs_c[d], descr_Ys_c[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'B');
            if(trsvtrsm == 'C' && spdnfactor == 'S' && factrs1 == 'U') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Us_sp1[d], descr_Xs_c[d], descr_Ys_c[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'B');
            if(trsvtrsm == 'R' && spdnfactor == 'S' && factrs1 == 'L') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Ls_sp1[d], descr_Xs_r[d], descr_Ys_r[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'B');
            if(trsvtrsm == 'R' && spdnfactor == 'S' && factrs1 == 'U') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Us_sp1[d], descr_Xs_r[d], descr_Ys_r[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'B');
            if(trsvtrsm == 'V' && spdnfactor == 'S' && factrs2syrk == 'L') gpu::spblas::trsv<T,I>(hs, 'H', descr_Ls_sp2[d], descr_Ys_vecs[d][0], descr_Xs_vecs[d][0], descrs_sparse_trsv2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'B');
            if(trsvtrsm == 'V' && spdnfactor == 'S' && factrs2syrk == 'U') gpu::spblas::trsv<T,I>(hs, 'N', descr_Us_sp2[d], descr_Ys_vecs[d][0], descr_Xs_vecs[d][0], descrs_sparse_trsv2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'B');
            if(trsvtrsm == 'C' && spdnfactor == 'S' && factrs2syrk == 'L') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Ls_sp2[d], descr_Ys_c[d], descr_Xs_c[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'B');
            if(trsvtrsm == 'C' && spdnfactor == 'S' && factrs2syrk == 'U') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Us_sp2[d], descr_Ys_c[d], descr_Xs_c[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'B');
            if(trsvtrsm == 'R' && spdnfactor == 'S' && factrs2syrk == 'L') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Ls_sp2[d], descr_Ys_r[d], descr_Xs_r[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'B');
            if(trsvtrsm == 'R' && spdnfactor == 'S' && factrs2syrk == 'U') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Us_sp2[d], descr_Ys_r[d], descr_Xs_r[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'B');
            if(factrs2syrk != 'S' && trsvtrsm != 'R') gpu::spblas::mm<T,I>(hs, 'N', 'N', descr_Bperms_sp[d], descr_Xs_c[d], descr_Fs_c[d], buffersizes_spmm[d], buffers_spmm[d], 'B');
            if(factrs2syrk != 'S' && trsvtrsm == 'R') gpu::spblas::mm<T,I>(hs, 'N', 'N', descr_Bperms_sp[d], descr_Xs_r[d], descr_Fs_r[d], buffersizes_spmm[d], buffers_spmm[d], 'B');

            gpu::dnblas::buffer_collect_size(hd, buffer_requirements.emplace_back(), [&](){
                T * dummyptrT = reinterpret_cast<T*>(sizeof(T));
                if(trsvtrsm == 'V' && spdnfactor == 'D' && factrs1 == 'L') gpu::dnblas::trsv<T,I>(hd, 'U', 'H', n_dofs_domain, ld_domain, dummyptrT, dummyptrT);
                if(trsvtrsm == 'V' && spdnfactor == 'D' && factrs1 == 'U') gpu::dnblas::trsv<T,I>(hd, 'L', 'N', n_dofs_domain, ld_domain, dummyptrT, dummyptrT);
                if(trsvtrsm == 'C' && spdnfactor == 'D' && factrs1 == 'L') gpu::dnblas::trsm<T,I>(hd, 'L', 'U', 'H', n_dofs_domain,    n_dofs_interface, dummyptrT, ld_domain, dummyptrT, ld_domain);
                if(trsvtrsm == 'C' && spdnfactor == 'D' && factrs1 == 'U') gpu::dnblas::trsm<T,I>(hd, 'L', 'L', 'N', n_dofs_domain,    n_dofs_interface, dummyptrT, ld_domain, dummyptrT, ld_domain);
                if(trsvtrsm == 'R' && spdnfactor == 'D' && factrs1 == 'L') gpu::dnblas::trsm<T,I>(hd, 'R', 'U', 'N', n_dofs_interface, n_dofs_domain,    dummyptrT, ld_domain, dummyptrT, ld_interface);
                if(trsvtrsm == 'R' && spdnfactor == 'D' && factrs1 == 'U') gpu::dnblas::trsm<T,I>(hd, 'R', 'L', 'H', n_dofs_interface, n_dofs_domain,    dummyptrT, ld_domain, dummyptrT, ld_interface);
                if(trsvtrsm == 'V' && spdnfactor == 'D' && factrs2syrk == 'L') gpu::dnblas::trsv<T,I>(hd, 'U', 'N', n_dofs_domain, ld_domain, dummyptrT, dummyptrT);
                if(trsvtrsm == 'V' && spdnfactor == 'D' && factrs2syrk == 'U') gpu::dnblas::trsv<T,I>(hd, 'L', 'H', n_dofs_domain, ld_domain, dummyptrT, dummyptrT);
                if(trsvtrsm == 'C' && spdnfactor == 'D' && factrs2syrk == 'L') gpu::dnblas::trsm<T,I>(hd, 'L', 'U', 'N', n_dofs_domain,    n_dofs_interface, dummyptrT, ld_domain, dummyptrT, ld_domain);
                if(trsvtrsm == 'C' && spdnfactor == 'D' && factrs2syrk == 'U') gpu::dnblas::trsm<T,I>(hd, 'L', 'L', 'H', n_dofs_domain,    n_dofs_interface, dummyptrT, ld_domain, dummyptrT, ld_domain);
                if(trsvtrsm == 'R' && spdnfactor == 'D' && factrs2syrk == 'L') gpu::dnblas::trsm<T,I>(hd, 'R', 'U', 'H', n_dofs_interface, n_dofs_domain,    dummyptrT, ld_domain, dummyptrT, ld_interface);
                if(trsvtrsm == 'R' && spdnfactor == 'D' && factrs2syrk == 'U') gpu::dnblas::trsm<T,I>(hd, 'R', 'L', 'N', n_dofs_interface, n_dofs_domain,    dummyptrT, ld_domain, dummyptrT, ld_interface);
                if(factrs2syrk == 'S' && trsvtrsm != 'R') gpu::dnblas::herk(hd, 'L', 'H', n_dofs_interface, n_dofs_domain, dummyptrT, ld_domain,    dummyptrT, ld_interface);
                if(factrs2syrk == 'S' && trsvtrsm == 'R') gpu::dnblas::herk(hd, 'L', 'N', n_dofs_interface, n_dofs_domain, dummyptrT, ld_interface, dummyptrT, ld_interface);
                gpu::dnblas::hemv(hd, 'L', n_dofs_interface, dummyptrT, ld_interface, dummyptrT, dummyptrT);
            });

            gpu::mgm::queue_wait(q);

            buffersizes_other[d] = *std::max_element(buffer_requirements.begin(), buffer_requirements.end());
        }
        tm_buffersize.stop();

        // permanent allocations for the lifetime of the object
        tm_alloc.start();
        {
            // host pinned memory
            tm_alloc_host.start();
            if(need_h_Us_sp)         h_Us_sp[d].resize(n_dofs_domain,    n_dofs_domain, n_nz_factor);
            if(need_h_Ls_sp)         h_Ls_sp[d].resize(n_dofs_domain,    n_dofs_domain, n_nz_factor);
            if(need_h_Bperms_sp) h_Bperms_sp[d].resize(n_dofs_interface, n_dofs_domain, feti.B1[d].nnz);
            if(applyalg == 'C')  h_applyc_xs[d].resize(n_dofs_interface);
            if(applyalg == 'C')  h_applyc_ys[d].resize(n_dofs_interface);
            tm_alloc_host.stop();

            // device memory
            tm_alloc_device.start();
            if(need_d_Us_sp)         d_Us_sp[d].resize(n_dofs_domain,    n_dofs_domain,    n_nz_factor);
            if(need_d_Ls_sp)         d_Ls_sp[d].resize(n_dofs_domain,    n_dofs_domain,    n_nz_factor);
            if(need_d_Bperms_sp) d_Bperms_sp[d].resize(n_dofs_interface, n_dofs_domain,    feti.B1[d].nnz);
            if(true)                    d_Fs[d].resize(n_dofs_interface, n_dofs_interface, ld_interface);
            d_apply_xs[d].resize(n_dofs_interface);
            d_apply_ys[d].resize(n_dofs_interface);
            d_applyg_D2Cs[d].resize(n_dofs_interface);
            if(spdnfactor == 'S') buffers_sptrs1[d] = gpu::mgm::memalloc_device(buffersizes_sptrs1[d]);
            if(spdnfactor == 'S' && factrs2syrk != 'S') buffers_sptrs2[d] = gpu::mgm::memalloc_device(buffersizes_sptrs2[d]);
            if(factrs2syrk != 'S') buffers_spmm[d] = gpu::mgm::memalloc_device(buffersizes_spmm[d]);
            tm_alloc_device.stop();
        }
        tm_alloc.stop();

        // set the pointers inside the descriptors of some matrices
        tm_setpointers.start();
        {
            if(need_d_Us_sp)     gpu::spblas::descr_matrix_csr_link_data(   descr_Us_sp1[d],     d_Us_sp[d]);
            if(need_d_Us_sp)     gpu::spblas::descr_matrix_csr_link_data(   descr_Us_sp2[d],     d_Us_sp[d]);
            if(need_d_Ls_sp)     gpu::spblas::descr_matrix_csr_link_data(   descr_Ls_sp1[d],     d_Ls_sp[d]);
            if(need_d_Ls_sp)     gpu::spblas::descr_matrix_csr_link_data(   descr_Ls_sp2[d],     d_Ls_sp[d]);
            if(need_d_Bperms_sp) gpu::spblas::descr_matrix_csr_link_data(descr_Bperms_sp[d], d_Bperms_sp[d]);
            gpu::spblas::descr_matrix_dense_link_data(descr_Fs_r[d], d_Fs[d]);
            gpu::spblas::descr_matrix_dense_link_data(descr_Fs_c[d], d_Fs[d]);
        }
        tm_setpointers.stop();

        // prepare matrices on host
        tm_Bperm.start();
        {
            Permutation<I> perm;
            perm.resize(n_dofs_domain);
            solvers_Kreg[d].getPermutation(perm);
            math::permuteColumns(h_Bperms_sp[d], feti.B1[d], perm);
            perm.clear();
        }
        tm_Bperm.stop();

        // extract symbolic pattern from the factor
        tm_get_factors.start();
        {
            tm_extract.start();
            Solver_Factors sym = DirectSparseSolver<T,I>::factorsSymmetry();
            if(sym == Solver_Factors::HERMITIAN_LOWER || sym == Solver_Factors::NONSYMMETRIC_BOTH) solvers_Kreg[d].getFactorL(h_Ls_sp[d], true, false);
            if(sym == Solver_Factors::HERMITIAN_UPPER || sym == Solver_Factors::NONSYMMETRIC_BOTH) solvers_Kreg[d].getFactorU(h_Us_sp[d], true, false);
            tm_extract.stop();
            tm_transpose.start();
            if(need_d_Ls_sp && sym == Solver_Factors::HERMITIAN_UPPER)
            {
                transmaps_U2L[d].resize(h_Us_sp[d].nnz);
                math::conjTransposeMapSetup(h_Ls_sp[d], transmaps_U2L[d], h_Us_sp[d]);
            }
            if(need_d_Us_sp && sym == Solver_Factors::HERMITIAN_LOWER)
            {
                transmaps_L2U[d].resize(h_Ls_sp[d].nnz);
                math::conjTransposeMapSetup(h_Us_sp[d], transmaps_L2U[d], h_Ls_sp[d]);
            }
            tm_transpose.stop();
        }
        tm_get_factors.stop();

        // copy some matrices to device
        tm_copyin.start();
        {
            if(need_d_Us_sp) gpu::mgm::copy_submit_h2d(q, d_Us_sp[d], h_Us_sp[d], true, false);
            if(need_d_Ls_sp) gpu::mgm::copy_submit_h2d(q, d_Ls_sp[d], h_Ls_sp[d], true, false);
            if(need_d_Bperms_sp) gpu::mgm::copy_submit_h2d(q, d_Bperms_sp[d], h_Bperms_sp[d], true, true);
            if(wcpset == 'W') gpu::mgm::queue_wait(q);
        }
        tm_copyin.stop();
        
        // proprocessing stage of the kernels
        tm_kernels_preprocess.start();
        {
            tm_trs1.start();
            if(trsvtrsm == 'V' && spdnfactor == 'S' && factrs1 == 'L') gpu::spblas::trsv<T,I>(hs, 'N', descr_Ls_sp1[d], descr_Xs_vecs[d][0], descr_Ys_vecs[d][0], descrs_sparse_trsv1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'P');
            if(trsvtrsm == 'V' && spdnfactor == 'S' && factrs1 == 'U') gpu::spblas::trsv<T,I>(hs, 'H', descr_Us_sp1[d], descr_Xs_vecs[d][0], descr_Ys_vecs[d][0], descrs_sparse_trsv1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'P');
            if(trsvtrsm == 'C' && spdnfactor == 'S' && factrs1 == 'L') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Ls_sp1[d], descr_Xs_c[d], descr_Ys_c[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'P');
            if(trsvtrsm == 'C' && spdnfactor == 'S' && factrs1 == 'U') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Us_sp1[d], descr_Xs_c[d], descr_Ys_c[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'P');
            if(trsvtrsm == 'R' && spdnfactor == 'S' && factrs1 == 'L') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Ls_sp1[d], descr_Xs_r[d], descr_Ys_r[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'P');
            if(trsvtrsm == 'R' && spdnfactor == 'S' && factrs1 == 'U') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Us_sp1[d], descr_Xs_r[d], descr_Ys_r[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'P');
            if(wcpset == 'W') gpu::mgm::queue_wait(q);
            tm_trs1.stop();

            if(factrs2syrk != 'S')
            {
                tm_trs2.start();
                if(trsvtrsm == 'V' && spdnfactor == 'S' && factrs2syrk == 'L') gpu::spblas::trsv<T,I>(hs, 'H', descr_Ls_sp2[d], descr_Ys_vecs[d][0], descr_Xs_vecs[d][0], descrs_sparse_trsv2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'P');
                if(trsvtrsm == 'V' && spdnfactor == 'S' && factrs2syrk == 'U') gpu::spblas::trsv<T,I>(hs, 'N', descr_Us_sp2[d], descr_Ys_vecs[d][0], descr_Xs_vecs[d][0], descrs_sparse_trsv2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'P');
                if(trsvtrsm == 'C' && spdnfactor == 'S' && factrs2syrk == 'L') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Ls_sp2[d], descr_Ys_c[d], descr_Xs_c[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'P');
                if(trsvtrsm == 'C' && spdnfactor == 'S' && factrs2syrk == 'U') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Us_sp2[d], descr_Ys_c[d], descr_Xs_c[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'P');
                if(trsvtrsm == 'R' && spdnfactor == 'S' && factrs2syrk == 'L') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Ls_sp2[d], descr_Ys_r[d], descr_Xs_r[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'P');
                if(trsvtrsm == 'R' && spdnfactor == 'S' && factrs2syrk == 'U') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Us_sp2[d], descr_Ys_r[d], descr_Xs_r[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'P');
                if(wcpset == 'W') gpu::mgm::queue_wait(q);
                tm_trs2.stop();
                
                tm_gemm.start();
                if(trsvtrsm != 'R') gpu::spblas::mm<T,I>(hs, 'N', 'N', descr_Bperms_sp[d], descr_Xs_c[d], descr_Fs_c[d], buffersizes_spmm[d], buffers_spmm[d], 'P');
                if(trsvtrsm == 'R') gpu::spblas::mm<T,I>(hs, 'N', 'N', descr_Bperms_sp[d], descr_Xs_r[d], descr_Fs_r[d], buffersizes_spmm[d], buffers_spmm[d], 'P');
                if(wcpset == 'W') gpu::mgm::queue_wait(q);
                tm_gemm.stop();
            }
        }
        tm_kernels_preprocess.stop();

        tm_mainloop_inner.stop();
    }
    tm_mainloop_outer.stop();

    // clean up the mess from buggy openmp in clang
    utils::run_dummy_parallel_region();
    
    // some stuff needed for apply
    tm_applystuff.start();
    d_applyg_x_cluster.resize(feti.lambdas.size);
    d_applyg_y_cluster.resize(feti.lambdas.size);
    {
        Vector_Dense<T*,I,Ah> h_applyg_xs_pointers;
        Vector_Dense<T*,I,Ah> h_applyg_ys_pointers;
        Vector_Dense<I,I,Ah> h_applyg_n_dofs_interfaces;
        Vector_Dense<I*,I,Ah> h_applyg_D2Cs_pointers;
        h_applyg_xs_pointers.resize(n_domains);
        h_applyg_ys_pointers.resize(n_domains);
        h_applyg_n_dofs_interfaces.resize(n_domains);
        h_applyg_D2Cs_pointers.resize(n_domains);
        d_applyg_xs_pointers.resize(n_domains);
        d_applyg_ys_pointers.resize(n_domains);
        d_applyg_n_dofs_interfaces.resize(n_domains);
        d_applyg_D2Cs_pointers.resize(n_domains);
        for(size_t d = 0; d < n_domains; d++)
        {
            h_applyg_xs_pointers.vals[d] = d_apply_xs[d].vals;
            h_applyg_ys_pointers.vals[d] = d_apply_ys[d].vals;
            h_applyg_n_dofs_interfaces.vals[d] = feti.B1[d].nrows;
            h_applyg_D2Cs_pointers.vals[d] = d_applyg_D2Cs[d].vals;
        }
        gpu::mgm::copy_submit_h2d(main_q, d_applyg_xs_pointers,       h_applyg_xs_pointers);
        gpu::mgm::copy_submit_h2d(main_q, d_applyg_ys_pointers,       h_applyg_ys_pointers);
        gpu::mgm::copy_submit_h2d(main_q, d_applyg_n_dofs_interfaces, h_applyg_n_dofs_interfaces);
        gpu::mgm::copy_submit_h2d(main_q, d_applyg_D2Cs_pointers,     h_applyg_D2Cs_pointers);
        for(size_t d = 0; d < n_domains; d++) gpu::mgm::copy_submit_h2d(main_q, d_applyg_D2Cs[d].vals, feti.D2C[d].data(), feti.B1[d].nrows);
        if(wcpset == 'W') gpu::mgm::queue_wait(main_q);
    }
    tm_applystuff.stop();

    // memory pool alloc
    tm_poolalloc.start();
    {
        // determine the mempool size such that no blocking would occur
        size_t mem_pool_size_request = 0;
        for(size_t qidx = 0; qidx < n_queues; qidx++)
        {
            size_t mem_pool_size_request_queue = 0;
            for(size_t d = qidx; d < n_domains; d += n_queues)
            {
                size_t mem_pool_size_request_domain = 0;
                I n_dofs_domain = feti.B1[d].ncols;
                I n_dofs_interface = feti.B1[d].nrows;
                I ld_domain = ((n_dofs_domain - 1) / align_elem + 1) * align_elem;
                I ld_interface = ((n_dofs_interface - 1) / align_elem + 1) * align_elem;
                if(need_d_Us_dn) mem_pool_size_request_domain += Matrix_Dense<T,I>::memoryRequirement(n_dofs_domain,    n_dofs_domain,    ld_domain);
                if(need_d_Ls_dn) mem_pool_size_request_domain += Matrix_Dense<T,I>::memoryRequirement(n_dofs_domain,    n_dofs_domain,    ld_domain);
                if(need_d_Xs_c)  mem_pool_size_request_domain += Matrix_Dense<T,I>::memoryRequirement(n_dofs_interface, n_dofs_domain,    ld_domain);
                if(need_d_Xs_r)  mem_pool_size_request_domain += Matrix_Dense<T,I>::memoryRequirement(n_dofs_domain,    n_dofs_interface, ld_interface);
                if(need_d_Ys_c)  mem_pool_size_request_domain += Matrix_Dense<T,I>::memoryRequirement(n_dofs_interface, n_dofs_domain,    ld_domain);
                if(need_d_Ys_r)  mem_pool_size_request_domain += Matrix_Dense<T,I>::memoryRequirement(n_dofs_domain,    n_dofs_interface, ld_interface);
                mem_pool_size_request_domain += buffersizes_other[d];
                // mem_pool_size_request_queue = std::max(mem_pool_size_request_queue, mem_pool_size_request_domain);
                mem_pool_size_request_queue += mem_pool_size_request_domain;
            }
            mem_pool_size_request += mem_pool_size_request_queue;
        }
        mem_pool_size_request = 2 * mem_pool_size_request + 1024; // safety margin
        size_t pool_size_device;
        gpu::mgm::memalloc_device_max(mem_pool_device, pool_size_device, mem_pool_size_request);
        // eslog::info("Mem pool requesting %zu MiB, allocated %zu MiB   (%zu B, %zu B)\n", mem_pool_size_request >> 20, pool_size_device >> 20, mem_pool_size_request, pool_size_device);
        size_t pool_size_device_aligned = (pool_size_device / align_B) * align_B;
        cbmba_res_device = std::make_unique<cbmba_resource>(mem_pool_device, pool_size_device_aligned);
    }
    tm_poolalloc.stop();

    tm_wait.start();
    gpu::mgm::device_wait();
    tm_wait.stop();
    tm_total.stop();

    stage = 2;

    print_timer("Set     total", tm_total);
    print_timer("Set       gpuinit", tm_gpuinit);
    print_timer("Set       gpuset", tm_gpuset);
    print_timer("Set       vecresize", tm_vecresize);
    print_timer("Set       gpucreate", tm_gpucreate);
    print_timer("Set       mainloop_outer", tm_mainloop_outer);
    print_timer("Set         mainloop_inner", tm_mainloop_inner);
    print_timer("Set           Kreg_combine", tm_Kreg_combine);
    print_timer("Set           solver_commit", tm_solver_commit);
    print_timer("Set           fact_symbolic", tm_fact_symbolic);
    print_timer("Set           descriptors", tm_descriptors);
    print_timer("Set           buffersize", tm_buffersize);
    print_timer("Set           alloc", tm_alloc);
    print_timer("Set             alloc_host", tm_alloc_host);
    print_timer("Set             alloc_device", tm_alloc_device);
    print_timer("Set           setpointers", tm_setpointers);
    print_timer("Set           Bperm", tm_Bperm);
    print_timer("Set           get_factors", tm_get_factors);
    print_timer("Set             extract", tm_extract);
    print_timer("Set             transpose", tm_transpose);
    print_timer("Set           copyin", tm_copyin);
    print_timer("Set           kernels_preprocess", tm_kernels_preprocess);
    print_timer("Set             trs1", tm_trs1);
    print_timer("Set             trs2", tm_trs2);
    print_timer("Set             gemm", tm_gemm);
    print_timer("Set       applystuff", tm_applystuff);
    print_timer("Set       poolalloc", tm_poolalloc);
    print_timer("Set       wait", tm_wait);
}

template <typename T, typename I>
void TotalFETIExplicitAcc<T,I>::update(const step::Step &step)
{
    if(stage != 2 && stage != 3) eslog::error("update: invalud order of operations in dualop\n");
    if(solvers_Kreg.size() != n_domains) eslog::error("non-matching number of matrices Kreg between set and update\n");
    if(DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::NONSYMMETRIC_BOTH && (factrs1 != 'L' || factrs2syrk != 'U')) eslog::error("Wrong options for non-symmetric K\n");

    my_timer tm_total, tm_mainloop_outer, tm_mainloop_inner, tm_Kreg_combine, tm_solver_commit, tm_fact_numeric, tm_get_factors, tm_extract, tm_transpose, tm_allocinpool, tm_setpointers, tm_copyin, tm_descr_update, tm_sp2dn, tm_kernels_compute, tm_trs1, tm_trs2, tm_gemm, tm_syrk, tm_freeinpool, tm_freeinpool_exec, tm_compute_d, tm_wait;

    tm_total.start();

    const bool need_L = (factrs1 == 'L' || factrs2syrk == 'L');
    const bool need_U = (factrs1 == 'U' || factrs2syrk == 'U');
    const bool need_d_Ls_dn = (need_L && spdnfactor == 'D');
    const bool need_d_Us_dn = (need_U && spdnfactor == 'D');
    const bool need_d_Ls_sp = need_L;
    const bool need_d_Us_sp = need_U;
    const bool need_d_Xs_c = (trsvtrsm == 'C' || trsvtrsm == 'V');
    const bool need_d_Xs_r = (trsvtrsm == 'R');
    const bool need_d_Ys_c = (need_d_Xs_c && spdnfactor == 'S');
    const bool need_d_Ys_r = (need_d_Xs_r && spdnfactor == 'S');
    const bool need_d_Xs_vecs = (trsvtrsm == 'V');
    const bool need_d_Ys_vecs = (need_d_Xs_vecs && spdnfactor == 'S');

    gpu::mgm::set_device(device);

    tm_mainloop_outer.start();
    #pragma omp parallel for schedule(static,1) if(wcpupdate == 'P')
    for(size_t d = 0; d < n_domains; d++)
    {
        tm_mainloop_inner.start();

        I n_dofs_domain = h_Bperms_sp[d].ncols;
        I n_dofs_interface = h_Bperms_sp[d].nrows;
        I ld_domain = ((n_dofs_domain - 1) / align_elem + 1) * align_elem;
        I ld_interface = ((n_dofs_interface - 1) / align_elem + 1) * align_elem;

        gpu::mgm::queue & q = queues[d % n_queues];
        gpu::dnblas::handle & hd = handles_dense[d % n_queues];
        gpu::spblas::handle & hs = handles_sparse[d % n_queues];

        void * buffer_other = nullptr;

        // Kreg = K + RegMat numeric values
        tm_Kreg_combine.start();
        {
            math::sumCombined(Kregs[d], T{1.0}, feti.K[d], feti.RegMat[d]);
        }
        tm_Kreg_combine.stop();

        // commit Kreg to solver (with numeric values)
        tm_solver_commit.start();
        {
            solvers_Kreg[d].commit(Kregs[d]);
        }
        tm_solver_commit.stop();
        
        // numeric factorization
        tm_fact_numeric.start();
        {
            solvers_Kreg[d].numericalFactorization();
        }
        tm_fact_numeric.stop();

        // extract values from numeric factor
        tm_get_factors.start();
        {
            tm_extract.start();
            Solver_Factors sym = DirectSparseSolver<T,I>::factorsSymmetry();
            if(sym == Solver_Factors::HERMITIAN_LOWER || sym == Solver_Factors::NONSYMMETRIC_BOTH) solvers_Kreg[d].getFactorL(h_Ls_sp[d], false, true);
            if(sym == Solver_Factors::HERMITIAN_UPPER || sym == Solver_Factors::NONSYMMETRIC_BOTH) solvers_Kreg[d].getFactorU(h_Us_sp[d], false, true);
            tm_extract.stop();
            tm_transpose.start();
            if(need_d_Ls_sp && sym == Solver_Factors::HERMITIAN_UPPER) math::conjTransposeMapUse(h_Ls_sp[d], transmaps_U2L[d], h_Us_sp[d]);
            if(need_d_Us_sp && sym == Solver_Factors::HERMITIAN_LOWER) math::conjTransposeMapUse(h_Us_sp[d], transmaps_L2U[d], h_Ls_sp[d]);
            tm_transpose.stop();
        }
        tm_get_factors.stop();

        // temporary allocations using the memory pool
        tm_allocinpool.start();
        {
            cbmba_d ator_dt(*cbmba_res_device, align_B);
            cbmba_res_device->do_transaction([&](){
                if(need_d_Us_dn) d_Us_dn[d] = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_dt);
                if(need_d_Ls_dn) d_Ls_dn[d] = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_dt);
                if(need_d_Xs_c)   d_Xs_c[d] = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_dt);
                if(need_d_Xs_r)   d_Xs_r[d] = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_dt);
                if(need_d_Ys_c)   d_Ys_c[d] = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_dt);
                if(need_d_Ys_r)   d_Ys_r[d] = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_dt);
                
                if(need_d_Us_dn) d_Us_dn[d]->resize(n_dofs_domain,    n_dofs_domain,    ld_domain);
                if(need_d_Ls_dn) d_Ls_dn[d]->resize(n_dofs_domain,    n_dofs_domain,    ld_domain);
                if(need_d_Xs_c)   d_Xs_c[d]->resize(n_dofs_interface, n_dofs_domain,    ld_domain);
                if(need_d_Xs_r)   d_Xs_r[d]->resize(n_dofs_domain,    n_dofs_interface, ld_interface);
                if(need_d_Ys_c)   d_Ys_c[d]->resize(n_dofs_interface, n_dofs_domain,    ld_domain);
                if(need_d_Ys_r)   d_Ys_r[d]->resize(n_dofs_domain,    n_dofs_interface, ld_interface);
                buffer_other = cbmba_res_device->allocate(buffersizes_other[d], align_B);
            });
        }
        tm_allocinpool.stop();

        if(spdnfactor == 'D') gpu::dnblas::buffer_set(hd, buffer_other, buffersizes_other[d]);

        // set the pointers inside the descriptors of the rest of the matrices
        tm_setpointers.start();
        {
            if(need_d_Us_dn) gpu::spblas::descr_matrix_dense_link_data(descr_Us_dn[d], *d_Us_dn[d]);
            if(need_d_Ls_dn) gpu::spblas::descr_matrix_dense_link_data(descr_Ls_dn[d], *d_Ls_dn[d]);
            if(need_d_Xs_c)  gpu::spblas::descr_matrix_dense_link_data(descr_Xs_c[d], *d_Xs_c[d]);
            if(need_d_Xs_r)  gpu::spblas::descr_matrix_dense_link_data(descr_Xs_r[d], *d_Xs_r[d]);
            if(need_d_Ys_c)  gpu::spblas::descr_matrix_dense_link_data(descr_Ys_c[d], *d_Ys_c[d]);
            if(need_d_Ys_r)  gpu::spblas::descr_matrix_dense_link_data(descr_Ys_r[d], *d_Ys_r[d]);
            if(need_d_Xs_vecs) for(I j = 0; j < n_dofs_interface; j++) gpu::spblas::descr_vector_dense_link_data(descr_Xs_vecs[d][j], *d_Xs_c[d], j);
            if(need_d_Ys_vecs) for(I j = 0; j < n_dofs_interface; j++) gpu::spblas::descr_vector_dense_link_data(descr_Ys_vecs[d][j], *d_Ys_c[d], j);
        }
        tm_setpointers.stop();

        // copy the new factors to device
        tm_copyin.start();
        {
            if(need_d_Us_sp) gpu::mgm::copy_submit_h2d(q, d_Us_sp[d], h_Us_sp[d], false, true);
            if(need_d_Ls_sp) gpu::mgm::copy_submit_h2d(q, d_Ls_sp[d], h_Ls_sp[d], false, true);
            if(wcpupdate == 'W') gpu::mgm::queue_wait(q);
        }
        tm_copyin.stop();

        // update sparse trsv/trsm descriptors to reflect the new matrix values
        tm_descr_update.start();
        {
            if(spdnfactor == 'S' && trsvtrsm == 'V' && factrs1 == 'L') gpu::spblas::trsv<T,I>(hs, 'N', descr_Ls_sp1[d], descr_Xs_vecs[d][0], descr_Ys_vecs[d][0], descrs_sparse_trsv1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'U');
            if(spdnfactor == 'S' && trsvtrsm == 'V' && factrs1 == 'U') gpu::spblas::trsv<T,I>(hs, 'H', descr_Us_sp1[d], descr_Xs_vecs[d][0], descr_Ys_vecs[d][0], descrs_sparse_trsv1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'U');
            if(spdnfactor == 'S' && trsvtrsm == 'C' && factrs1 == 'L') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Ls_sp1[d], descr_Xs_c[d], descr_Ys_c[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'U');
            if(spdnfactor == 'S' && trsvtrsm == 'C' && factrs1 == 'U') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Us_sp1[d], descr_Xs_c[d], descr_Ys_c[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'U');
            if(spdnfactor == 'S' && trsvtrsm == 'R' && factrs1 == 'L') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Ls_sp1[d], descr_Xs_r[d], descr_Ys_r[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'U');
            if(spdnfactor == 'S' && trsvtrsm == 'R' && factrs1 == 'U') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Us_sp1[d], descr_Xs_r[d], descr_Ys_r[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'U');
            if(spdnfactor == 'S' && trsvtrsm == 'V' && factrs2syrk == 'L') gpu::spblas::trsv<T,I>(hs, 'H', descr_Ls_sp2[d], descr_Ys_vecs[d][0], descr_Xs_vecs[d][0], descrs_sparse_trsv2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'U');
            if(spdnfactor == 'S' && trsvtrsm == 'V' && factrs2syrk == 'U') gpu::spblas::trsv<T,I>(hs, 'N', descr_Us_sp2[d], descr_Ys_vecs[d][0], descr_Xs_vecs[d][0], descrs_sparse_trsv2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'U');
            if(spdnfactor == 'S' && trsvtrsm == 'C' && factrs2syrk == 'L') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Ls_sp2[d], descr_Ys_c[d], descr_Xs_c[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'U');
            if(spdnfactor == 'S' && trsvtrsm == 'C' && factrs2syrk == 'U') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Us_sp2[d], descr_Ys_c[d], descr_Xs_c[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'U');
            if(spdnfactor == 'S' && trsvtrsm == 'R' && factrs2syrk == 'L') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Ls_sp2[d], descr_Ys_r[d], descr_Xs_r[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'U');
            if(spdnfactor == 'S' && trsvtrsm == 'R' && factrs2syrk == 'U') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Us_sp2[d], descr_Ys_r[d], descr_Xs_r[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'U');
            if(wcpupdate == 'W') gpu::mgm::queue_wait(q);
        }
        tm_descr_update.stop();
        
        // sparse to dense on device
        tm_sp2dn.start();
        {
            if(need_d_Us_dn) gpu::spblas::sparse_to_dense<T,I>(hs, 'N', descr_Us_sp1[d],    descr_Us_dn[d], buffersizes_other[d], buffer_other, 'C');
            if(need_d_Ls_dn) gpu::spblas::sparse_to_dense<T,I>(hs, 'N', descr_Ls_sp1[d],    descr_Ls_dn[d], buffersizes_other[d], buffer_other, 'C');
            if(need_d_Xs_c)  gpu::spblas::sparse_to_dense<T,I>(hs, 'T', descr_Bperms_sp[d], descr_Xs_c[d], buffersizes_other[d], buffer_other, 'C');
            if(need_d_Xs_r)  gpu::spblas::sparse_to_dense<T,I>(hs, 'T', descr_Bperms_sp[d], descr_Xs_r[d], buffersizes_other[d], buffer_other, 'C');
            if(wcpupdate == 'W') gpu::mgm::queue_wait(q);
        }
        tm_sp2dn.stop();

        // perform the actual assembly
        tm_kernels_compute.start();
        {
            tm_trs1.start();
            if(trsvtrsm == 'V' && spdnfactor == 'S' && factrs1 == 'L') for(I j = 0; j < n_dofs_interface; j++) gpu::spblas::trsv<T,I>(hs, 'N', descr_Ls_sp1[d], descr_Xs_vecs[d][j], descr_Ys_vecs[d][j], descrs_sparse_trsv1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'C');
            if(trsvtrsm == 'V' && spdnfactor == 'S' && factrs1 == 'U') for(I j = 0; j < n_dofs_interface; j++) gpu::spblas::trsv<T,I>(hs, 'H', descr_Us_sp1[d], descr_Xs_vecs[d][j], descr_Ys_vecs[d][j], descrs_sparse_trsv1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'C');
            if(trsvtrsm == 'C' && spdnfactor == 'S' && factrs1 == 'L') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Ls_sp1[d], descr_Xs_c[d], descr_Ys_c[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'C');
            if(trsvtrsm == 'C' && spdnfactor == 'S' && factrs1 == 'U') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Us_sp1[d], descr_Xs_c[d], descr_Ys_c[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'C');
            if(trsvtrsm == 'R' && spdnfactor == 'S' && factrs1 == 'L') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Ls_sp1[d], descr_Xs_r[d], descr_Ys_r[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'C');
            if(trsvtrsm == 'R' && spdnfactor == 'S' && factrs1 == 'U') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Us_sp1[d], descr_Xs_r[d], descr_Ys_r[d], descrs_sparse_trsm1[d], buffersizes_sptrs1[d], buffers_sptrs1[d], 'C');
            if(trsvtrsm == 'V' && spdnfactor == 'D' && factrs1 == 'L') for(I j = 0; j < n_dofs_interface; j++) gpu::dnblas::trsv<T,I>(hd, 'U', 'H', n_dofs_domain, ld_domain, d_Ls_dn[d]->vals, d_Xs_c[d]->vals + j * d_Xs_c[d]->get_ld());
            if(trsvtrsm == 'V' && spdnfactor == 'D' && factrs1 == 'U') for(I j = 0; j < n_dofs_interface; j++) gpu::dnblas::trsv<T,I>(hd, 'L', 'N', n_dofs_domain, ld_domain, d_Us_dn[d]->vals, d_Xs_c[d]->vals + j * d_Xs_c[d]->get_ld());
            if(trsvtrsm == 'C' && spdnfactor == 'D' && factrs1 == 'L') gpu::dnblas::trsm<T,I>(hd, 'L', 'U', 'H', n_dofs_domain,    n_dofs_interface, d_Ls_dn[d]->vals, ld_domain, d_Xs_c[d]->vals, ld_domain);
            if(trsvtrsm == 'C' && spdnfactor == 'D' && factrs1 == 'U') gpu::dnblas::trsm<T,I>(hd, 'L', 'L', 'N', n_dofs_domain,    n_dofs_interface, d_Us_dn[d]->vals, ld_domain, d_Xs_c[d]->vals, ld_domain);
            if(trsvtrsm == 'R' && spdnfactor == 'D' && factrs1 == 'L') gpu::dnblas::trsm<T,I>(hd, 'R', 'U', 'N', n_dofs_interface, n_dofs_domain,    d_Ls_dn[d]->vals, ld_domain, d_Xs_r[d]->vals, ld_interface);
            if(trsvtrsm == 'R' && spdnfactor == 'D' && factrs1 == 'U') gpu::dnblas::trsm<T,I>(hd, 'R', 'L', 'H', n_dofs_interface, n_dofs_domain,    d_Us_dn[d]->vals, ld_domain, d_Xs_r[d]->vals, ld_interface);
            if(wcpupdate == 'W') gpu::mgm::queue_wait(q);
            tm_trs1.stop();

            if(factrs2syrk == 'S')
            {
                tm_syrk.start();
                if(spdnfactor == 'S' && trsvtrsm != 'R') gpu::dnblas::herk(hd, 'L', 'H', n_dofs_interface, n_dofs_domain, d_Ys_c[d]->vals, ld_domain,    d_Fs[d].vals, ld_interface);
                if(spdnfactor == 'S' && trsvtrsm == 'R') gpu::dnblas::herk(hd, 'L', 'N', n_dofs_interface, n_dofs_domain, d_Ys_r[d]->vals, ld_interface, d_Fs[d].vals, ld_interface);
                if(spdnfactor == 'D' && trsvtrsm != 'R') gpu::dnblas::herk(hd, 'L', 'H', n_dofs_interface, n_dofs_domain, d_Xs_c[d]->vals, ld_domain,    d_Fs[d].vals, ld_interface);
                if(spdnfactor == 'D' && trsvtrsm == 'R') gpu::dnblas::herk(hd, 'L', 'N', n_dofs_interface, n_dofs_domain, d_Xs_r[d]->vals, ld_interface, d_Fs[d].vals, ld_interface);
                if(wcpupdate == 'W') gpu::mgm::queue_wait(q);
                tm_syrk.stop();
            }
            else
            {
                tm_trs2.start();
                if(trsvtrsm == 'V' && spdnfactor == 'S' && factrs2syrk == 'L') for(I j = 0; j < n_dofs_interface; j++) gpu::spblas::trsv<T,I>(hs, 'H', descr_Ls_sp2[d], descr_Ys_vecs[d][j], descr_Xs_vecs[d][j], descrs_sparse_trsv2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'C');
                if(trsvtrsm == 'V' && spdnfactor == 'S' && factrs2syrk == 'U') for(I j = 0; j < n_dofs_interface; j++) gpu::spblas::trsv<T,I>(hs, 'N', descr_Us_sp2[d], descr_Ys_vecs[d][j], descr_Xs_vecs[d][j], descrs_sparse_trsv2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'C');
                if(trsvtrsm == 'C' && spdnfactor == 'S' && factrs2syrk == 'L') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Ls_sp2[d], descr_Ys_c[d], descr_Xs_c[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'C');
                if(trsvtrsm == 'C' && spdnfactor == 'S' && factrs2syrk == 'U') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Us_sp2[d], descr_Ys_c[d], descr_Xs_c[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'C');
                if(trsvtrsm == 'R' && spdnfactor == 'S' && factrs2syrk == 'L') gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', descr_Ls_sp2[d], descr_Ys_r[d], descr_Xs_r[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'C');
                if(trsvtrsm == 'R' && spdnfactor == 'S' && factrs2syrk == 'U') gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', descr_Us_sp2[d], descr_Ys_r[d], descr_Xs_r[d], descrs_sparse_trsm2[d], buffersizes_sptrs2[d], buffers_sptrs2[d], 'C');
                if(trsvtrsm == 'V' && spdnfactor == 'D' && factrs2syrk == 'L') for(I j = 0; j < n_dofs_interface; j++) gpu::dnblas::trsv<T,I>(hd, 'U', 'N', n_dofs_domain, ld_domain, d_Ls_dn[d]->vals, d_Xs_c[d]->vals + j * d_Xs_c[d]->get_ld());
                if(trsvtrsm == 'V' && spdnfactor == 'D' && factrs2syrk == 'U') for(I j = 0; j < n_dofs_interface; j++) gpu::dnblas::trsv<T,I>(hd, 'L', 'H', n_dofs_domain, ld_domain, d_Us_dn[d]->vals, d_Xs_c[d]->vals + j * d_Xs_c[d]->get_ld());
                if(trsvtrsm == 'C' && spdnfactor == 'D' && factrs2syrk == 'L') gpu::dnblas::trsm<T,I>(hd, 'L', 'U', 'N', n_dofs_domain,    n_dofs_interface, d_Ls_dn[d]->vals, ld_domain, d_Xs_c[d]->vals, ld_domain);
                if(trsvtrsm == 'C' && spdnfactor == 'D' && factrs2syrk == 'U') gpu::dnblas::trsm<T,I>(hd, 'L', 'L', 'H', n_dofs_domain,    n_dofs_interface, d_Us_dn[d]->vals, ld_domain, d_Xs_c[d]->vals, ld_domain);
                if(trsvtrsm == 'R' && spdnfactor == 'D' && factrs2syrk == 'L') gpu::dnblas::trsm<T,I>(hd, 'R', 'U', 'H', n_dofs_interface, n_dofs_domain,    d_Ls_dn[d]->vals, ld_domain, d_Xs_r[d]->vals, ld_interface);
                if(trsvtrsm == 'R' && spdnfactor == 'D' && factrs2syrk == 'U') gpu::dnblas::trsm<T,I>(hd, 'R', 'L', 'N', n_dofs_interface, n_dofs_domain,    d_Us_dn[d]->vals, ld_domain, d_Xs_r[d]->vals, ld_interface);
                if(wcpupdate == 'W') gpu::mgm::queue_wait(q);
                tm_trs2.stop();

                tm_gemm.start();
                if(trsvtrsm != 'R') gpu::spblas::mm<T,I>(hs, 'N', 'N', descr_Bperms_sp[d], descr_Xs_c[d], descr_Fs_c[d], buffersizes_spmm[d], buffers_spmm[d], 'C');
                if(trsvtrsm == 'R') gpu::spblas::mm<T,I>(hs, 'N', 'N', descr_Bperms_sp[d], descr_Xs_r[d], descr_Fs_r[d], buffersizes_spmm[d], buffers_spmm[d], 'C');
                if(wcpupdate == 'W') gpu::mgm::queue_wait(q);
                tm_gemm.stop();
            }
        }
        tm_kernels_compute.stop();

        if(spdnfactor == 'D') gpu::dnblas::buffer_unset(hd);

        // free the temporary memory from the pool
        tm_freeinpool.start();
        gpu::mgm::submit_host_function(q, [&,d,buffer_other](){
            tm_freeinpool_exec.start();
            if(need_d_Us_dn) { d_Us_dn[d].reset(nullptr); }
            if(need_d_Ls_dn) { d_Ls_dn[d].reset(nullptr); }
            if(need_d_Xs_c)  { d_Xs_c[d].reset(nullptr); }
            if(need_d_Xs_r)  { d_Xs_r[d].reset(nullptr); }
            if(need_d_Ys_c)  { d_Ys_c[d].reset(nullptr); }
            if(need_d_Ys_r)  { d_Ys_r[d].reset(nullptr); }
            cbmba_res_device->deallocate(buffer_other);
            tm_freeinpool_exec.stop();
        });
        if(wcpupdate == 'W') gpu::mgm::queue_wait(q);
        tm_freeinpool.stop();

        tm_mainloop_inner.stop();
    }
    tm_mainloop_outer.stop();

    // clean up the mess from buggy openmp in clang
    utils::run_dummy_parallel_region();

    tm_compute_d.start();
    {
        // just use the cpu solver
        std::vector<Vector_Dense<T,I>> Kplus_fs(n_domains);
        #pragma omp parallel for schedule(static,1) if(wcpupdate == 'P')
        for(size_t d = 0; d < n_domains; d++) {
            solvers_Kreg[d].solve(feti.f[d], Kplus_fs[d]);
        }
        applyB(feti, Kplus_fs, d);
        math::add(d, T{-1}, feti.c);
    }
    tm_compute_d.stop();

    tm_wait.start();
    gpu::mgm::device_wait();
    tm_wait.stop();
    tm_total.stop();

    stage = 3;

    print_timer("Update  total", tm_total);
    print_timer("Update    mainloop_outer", tm_mainloop_outer);
    print_timer("Update      mainloop_inner", tm_mainloop_inner);
    print_timer("Update        Kreg_combine", tm_Kreg_combine);
    print_timer("Update        solver_commit", tm_solver_commit);
    print_timer("Update        fact_numeric", tm_fact_numeric);
    print_timer("Update        get_factors", tm_get_factors);
    print_timer("Update          extract", tm_extract);
    print_timer("Update          transpose", tm_transpose);
    print_timer("Update        allocinpool", tm_allocinpool);
    print_timer("Update        setpointers", tm_setpointers);
    print_timer("Update        copyin", tm_copyin);
    print_timer("Update        descr_update", tm_descr_update);
    print_timer("Update        sp2dn", tm_sp2dn);
    print_timer("Update        kernels_compute", tm_kernels_compute);
    print_timer("Update          trs1", tm_trs1);
    print_timer("Update          trs2", tm_trs2);
    print_timer("Update          gemm", tm_gemm);
    print_timer("Update          syrk", tm_syrk);
    print_timer("Update        freeinpool", tm_freeinpool);
    print_timer("Update          freeinpool_exec", tm_freeinpool_exec);
    print_timer("Update    compute_d", tm_compute_d);
    print_timer("Update    wait", tm_wait);

//    print(step);
}

template <typename T, typename I>
void TotalFETIExplicitAcc<T,I>::apply(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    if(stage != 3) eslog::error("invalid stage when calling apply\n");

    gpu::mgm::set_device(device);

    if(applyalg == 'G')
    {
        my_timer tm_total, tm_copyin, tm_scatter, tm_mv_outer, tm_mv, tm_zerofill, tm_gather, tm_copyout, tm_wait;

        tm_total.start();

        // copy x_cluster to device
        tm_copyin.start();
        gpu::mgm::copy_submit_h2d(main_q, d_applyg_x_cluster, x_cluster);
        if(wcpapply == 'W') gpu::mgm::queue_wait(main_q);
        tm_copyin.stop();

        // scatter
        tm_scatter.start();
        gpu::kernels::DCmap_scatter(main_q, d_applyg_xs_pointers, d_applyg_n_dofs_interfaces, d_applyg_x_cluster, d_applyg_D2Cs_pointers);
        if(wcpapply == 'W') gpu::mgm::queue_wait(main_q);
        tm_scatter.stop();

        // apply
        tm_mv_outer.start();
        gpu::mgm::queue_async_barrier({main_q}, queues);
        #pragma omp parallel for schedule(static,1) if(wcpapply == 'P')
        for(size_t d = 0; d < n_domains; d++)
        {
            gpu::mgm::queue & q = queues[d % n_queues];
            gpu::dnblas::handle & hd = handles_dense[d % n_queues];

            tm_mv.start();
            gpu::dnblas::hemv(hd, 'L', d_Fs[d].nrows, d_Fs[d].vals, d_Fs[d].get_ld(), d_apply_xs[d].vals, d_apply_ys[d].vals);
            if(wcpapply == 'W') gpu::mgm::queue_wait(q);
            tm_mv.stop();
        }
        tm_mv_outer.stop();

        // zerofill y_cluster on device
        tm_zerofill.start();
        gpu::mgm::memset_submit(main_q, d_applyg_y_cluster.vals, d_applyg_y_cluster.size * sizeof(T), 0);
        if(wcpapply == 'W') gpu::mgm::queue_wait(main_q);
        tm_zerofill.stop();

        gpu::mgm::queue_async_barrier(queues, {main_q});

        // gather
        tm_gather.start();
        gpu::kernels::DCmap_gather(main_q, d_applyg_ys_pointers, d_applyg_n_dofs_interfaces, d_applyg_y_cluster, d_applyg_D2Cs_pointers);
        if(wcpapply == 'W') gpu::mgm::queue_wait(main_q);
        tm_gather.stop();

        // copy y_cluster from device
        tm_copyout.start();
        gpu::mgm::copy_submit_d2h(main_q, y_cluster, d_applyg_y_cluster);
        if(wcpapply == 'W') gpu::mgm::queue_wait(main_q);
        tm_copyout.stop();

        // wait
        tm_wait.start();
        gpu::mgm::device_wait();
        tm_wait.stop();

        tm_total.stop();

        print_timer("Apply   total", tm_total);
        print_timer("Apply     copyin", tm_copyin);
        print_timer("Apply     scatter", tm_scatter);
        print_timer("Apply     mv_outer", tm_mv_outer);
        print_timer("Apply       mv", tm_mv);
        print_timer("Apply     zerofill", tm_zerofill);
        print_timer("Apply     gather", tm_gather);
        print_timer("Apply     copyout", tm_copyout);
        print_timer("Apply     wait", tm_wait);
    }
    
    if(applyalg == 'C')
    {
        my_timer tm_total, tm_apply_outer, tm_apply_inner, tm_scatter, tm_copyin, tm_mv, tm_copyout, tm_zerofill, tm_wait, tm_gather, tm_gather_inner;

        tm_total.start();

        tm_apply_outer.start();
        #pragma omp parallel for schedule(static,1) if(wcpapply == 'P')
        for(size_t d = 0; d < n_domains; d++)
        {
            tm_apply_inner.start();

            I n_dofs_interface = h_Bperms_sp[d].nrows;

            gpu::mgm::queue & q = queues[d % n_queues];
            gpu::dnblas::handle & hd = handles_dense[d % n_queues];

            tm_scatter.start();
            for(I i = 0; i < n_dofs_interface; i++) h_applyc_xs[d].vals[i] = x_cluster.vals[feti.D2C[d][i]];
            tm_scatter.stop();

            tm_copyin.start();
            gpu::mgm::copy_submit_h2d(q, d_apply_xs[d], h_applyc_xs[d]);
            if(wcpapply == 'W') gpu::mgm::queue_wait(q);
            tm_copyin.stop();

            tm_mv.start();
            gpu::dnblas::hemv(hd, 'L', d_Fs[d].nrows, d_Fs[d].vals, d_Fs[d].get_ld(), d_apply_xs[d].vals, d_apply_ys[d].vals);
            if(wcpapply == 'W') gpu::mgm::queue_wait(q);
            tm_mv.stop();

            tm_copyout.start();
            gpu::mgm::copy_submit_d2h(q, h_applyc_ys[d], d_apply_ys[d]);
            if(wcpapply == 'W') gpu::mgm::queue_wait(q);
            tm_copyout.stop();

            tm_apply_inner.stop();
        }
        tm_apply_outer.stop();

        tm_zerofill.start();
        std::fill_n(y_cluster.vals, y_cluster.Vector_Dense<T,I>::size, T{0});
        tm_zerofill.stop();

        tm_wait.start();
        gpu::mgm::device_wait();
        tm_wait.stop();

        tm_gather.start();
        #pragma omp parallel for schedule(static,1) if(wcpapply == 'P')
        for(size_t d = 0; d < n_domains; d++)
        {
            tm_gather_inner.start();
            I n_dofs_interface = h_Bperms_sp[d].nrows;
            for(I i = 0; i < n_dofs_interface; i++)
            {
                if constexpr(utils::is_real<T>())
                {
                    #pragma omp atomic
                    y_cluster.vals[feti.D2C[d][i]] += h_applyc_ys[d].vals[i];
                }
                if constexpr(utils::is_complex<T>())
                {
                    #pragma omp atomic
                    utils::real_ref(y_cluster.vals[feti.D2C[d][i]]) += utils::real_ref(h_applyc_ys[d].vals[i]);
                    #pragma omp atomic
                    utils::imag_ref(y_cluster.vals[feti.D2C[d][i]]) += utils::imag_ref(h_applyc_ys[d].vals[i]);
                }
            }
            tm_gather_inner.stop();
        }
        tm_gather.stop();

        tm_total.stop();

        print_timer("Apply   total", tm_total);
        print_timer("Apply     apply_outer", tm_apply_outer);
        print_timer("Apply       apply_inner", tm_apply_inner);
        print_timer("Apply         scatter", tm_scatter);
        print_timer("Apply         copyin", tm_copyin);
        print_timer("Apply         mv", tm_mv);
        print_timer("Apply         copyout", tm_copyout);
        print_timer("Apply     zerofill", tm_zerofill);
        print_timer("Apply     wait", tm_wait);
        print_timer("Apply     gather", tm_gather);
        print_timer("Apply     gather_inner", tm_gather_inner);
    }

    y_cluster.synchronize();
}

template <typename T, typename I>
void TotalFETIExplicitAcc<T,I>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    // just do it on cpu
    #pragma omp parallel for schedule(static,1) if(wcpset == 'P')
    for (size_t d = 0; d < feti.K.size(); ++d) {
        Vector_Dense<T,I> z;
        z.resize(y[d]);
        applyBt(feti, d, x, z, T{-1});
        math::add(z, T{1}, feti.f[d]);
        solvers_Kreg[d].solve(z, y[d]);
    }
}

template <typename T, typename I>
void TotalFETIExplicitAcc<T,I>::print(const step::Step &step)
{
    eslog::error("todo\n");
    // if (info::ecf->output.print_matrices) {
    //     eslog::storedata(" STORE: feti/dualop/{Kplus}\n");
    //     for (size_t di = 0; di < feti.K.size(); ++di) {
    //         math::store(Kplus[di], utils::filename(utils::debugDirectory(step) + "/feti/dualop", (std::string("Kplus") + std::to_string(di)).c_str()).c_str());
    //     }
    //     math::store(d, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "d").c_str());
    // }
}

//template class TotalFETIExplicitAcc<float , int32_t>;
template class TotalFETIExplicitAcc<double, int32_t>;
// template class TotalFETIExplicitAcc<float , int64_t>;
// template class TotalFETIExplicitAcc<double, int64_t>;
//template class TotalFETIExplicitAcc<std::complex<float >, int32_t>;
//template class TotalFETIExplicitAcc<std::complex<double>, int32_t>;
// template class TotalFETIExplicitAcc<std::complex<float >, int64_t>;
// template class TotalFETIExplicitAcc<std::complex<double>, int64_t>;

}
