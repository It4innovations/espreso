
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

using CONCURRENCY             = DualOperatorExplicitGpuConfig::CONCURRENCY;
using MATRIX_STORAGE          = DualOperatorExplicitGpuConfig::MATRIX_STORAGE;
using TRSM1_SOLVE_TYPE        = DualOperatorExplicitGpuConfig::TRSM1_SOLVE_TYPE;
using TRSM2_SOLVE_TYPE        = DualOperatorExplicitGpuConfig::TRSM2_SOLVE_TYPE;
using MATRIX_ORDER            = DualOperatorExplicitGpuConfig::MATRIX_ORDER;
using PATH_IF_HERMITIAN       = DualOperatorExplicitGpuConfig::PATH_IF_HERMITIAN;
using TRIANGLE_MATRIX_SHARING = DualOperatorExplicitGpuConfig::TRIANGLE_MATRIX_SHARING;
using QUEUE_COUNT             = DualOperatorExplicitGpuConfig::QUEUE_COUNT;
using DEVICE                  = DualOperatorExplicitGpuConfig::DEVICE;

template <typename T, typename I>
TotalFETIExplicitAcc<T,I>::TotalFETIExplicitAcc(FETI<T> &feti)
: DualOperator<T>(feti), n_domains(0), n_queues(0), mem_pool_device(nullptr)
{
    if(stage != 0) eslog::error("init: invalid order of operations in dualop\n");

    config = &feti.configuration.dual_operator_explicit_gpu_config;
    config_replace_defaults();
    if(config->trsm_rhs_sol_order == MATRIX_ORDER::ROW_MAJOR) order_X = 'R';
    if(config->trsm_rhs_sol_order == MATRIX_ORDER::COL_MAJOR) order_X = 'C';
    order_F = order_X;
    is_system_hermitian = (DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_LOWER || DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_UPPER);
    is_factor1_dense = (config->trsm1_factor_storage == MATRIX_STORAGE::DENSE);
    is_factor2_dense = (config->trsm2_factor_storage == MATRIX_STORAGE::DENSE);
    is_factor1_sparse = (config->trsm1_factor_storage == MATRIX_STORAGE::SPARSE);
    is_factor2_sparse = (config->trsm2_factor_storage == MATRIX_STORAGE::SPARSE);
    is_path_trsm = (config->path_if_hermitian == PATH_IF_HERMITIAN::TRSM);
    is_path_herk = (config->path_if_hermitian == PATH_IF_HERMITIAN::HERK);
    trsm1_use_L  = (config->trsm1_solve_type == TRSM1_SOLVE_TYPE::L);
    trsm1_use_LH = (config->trsm1_solve_type == TRSM1_SOLVE_TYPE::LHH);
    trsm2_use_U  = (config->trsm2_solve_type == TRSM2_SOLVE_TYPE::U && is_path_trsm);
    trsm2_use_UH = (config->trsm2_solve_type == TRSM2_SOLVE_TYPE::UHH && is_path_trsm);
    need_Y = (is_factor1_sparse || is_factor2_sparse);
    Solver_Factors sym = DirectSparseSolver<T,I>::factorsSymmetry();
    solver_get_L = (sym == Solver_Factors::HERMITIAN_LOWER || sym == Solver_Factors::NONSYMMETRIC_BOTH);
    solver_get_U = (sym == Solver_Factors::HERMITIAN_UPPER || sym == Solver_Factors::NONSYMMETRIC_BOTH);
    can_use_LH_is_U_h_sp = (is_system_hermitian && (trsm2_use_U || solver_get_U));
    can_use_UH_is_L_h_sp = (is_system_hermitian && (trsm1_use_L || solver_get_L));
    can_use_LH_is_U_d_sp = (is_system_hermitian &&  trsm2_use_U);
    can_use_UH_is_L_d_sp = (is_system_hermitian &&  trsm1_use_L);
    can_use_LH_is_U_d_dn = (is_system_hermitian &&  trsm2_use_U && is_factor2_dense);
    can_use_UH_is_L_d_dn = (is_system_hermitian &&  trsm1_use_L && is_factor1_dense);
    need_conjtrans_L2LH = (is_system_hermitian && solver_get_L && (trsm1_use_LH || trsm2_use_U)) || (!is_system_hermitian && trsm1_use_LH);
    need_conjtrans_U2UH = (is_system_hermitian && solver_get_U && (trsm2_use_UH || trsm1_use_L)) || (!is_system_hermitian && trsm2_use_UH);
    is_f_triangles_shared = (is_system_hermitian && config->f_sharing_if_hermitian == TRIANGLE_MATRIX_SHARING::SHARED);
    need_f_tmp = (is_f_triangles_shared && is_path_trsm);

    if(is_path_herk && !is_system_hermitian) eslog::error("cannot do herk path with non-hermitian system\n");

    device = gpu::mgm::get_device_by_mpi(info::mpi::rank, info::mpi::size);

    stage = 1;
}

template <typename T, typename I>
TotalFETIExplicitAcc<T,I>::~TotalFETIExplicitAcc()
{
    gpu::mgm::set_device(device);

    my_timer tm_total, tm_descriptors, tm_gpulibs, tm_gpumem, tm_gpumempool, tm_gpuhostmem, tm_wait;

    tm_total.start();

    tm_descriptors.start();
    for(size_t d = 0; d < n_domains; d++)
    {
        per_domain_stuff & data = domain_data[d];

        gpu::spblas::descr_matrix_csr_destroy(data.descr_L_sp);
        gpu::spblas::descr_matrix_csr_destroy(data.descr_LH_sp);
        gpu::spblas::descr_matrix_csr_destroy(data.descr_U_sp);
        gpu::spblas::descr_matrix_csr_destroy(data.descr_UH_sp);
        gpu::spblas::descr_matrix_csr_destroy(data.descr_Bperm_sp);
        gpu::spblas::descr_matrix_dense_destroy(data.descr_L_dn);
        gpu::spblas::descr_matrix_dense_destroy(data.descr_LH_dn);
        gpu::spblas::descr_matrix_dense_destroy(data.descr_U_dn);
        gpu::spblas::descr_matrix_dense_destroy(data.descr_UH_dn);
        gpu::spblas::descr_matrix_dense_destroy(data.descr_X);
        gpu::spblas::descr_matrix_dense_destroy(data.descr_Y);
        gpu::spblas::descr_matrix_dense_destroy(data.descr_F);
        gpu::spblas::descr_sparse_trsm_destroy(data.descr_sparse_trsm1);
        gpu::spblas::descr_sparse_trsm_destroy(data.descr_sparse_trsm2);
    }
    tm_descriptors.stop();

    tm_gpulibs.start();
    gpu::mgm::queue_destroy(main_q);
    for(gpu::mgm::queue & q : queues) gpu::mgm::queue_destroy(q);
    for(gpu::dnblas::handle & h : handles_dense) gpu::dnblas::handle_destroy(h);
    for(gpu::spblas::handle & h : handles_sparse) gpu::spblas::handle_destroy(h);
    tm_gpulibs.stop();

    tm_gpumem.start();
    tm_gpumempool.start();
    gpu::mgm::memfree_device(mem_pool_device);
    tm_gpumempool.stop();
    for(auto & f : d_Fs_allocated) f.clear();
    for(size_t d = 0; d < n_domains; d++)
    {
        per_domain_stuff & data = domain_data[d];
        data.d_L_sp.clear();
        data.d_LH_sp.clear();
        data.d_U_sp.clear();
        data.d_UH_sp.clear();
        data.d_Bperm_sp.clear();
        data.d_applyg_D2C.clear();
        data.d_apply_x.clear();
        data.d_apply_y.clear();
        gpu::mgm::memfree_device(data.buffer_sptrs1);
        gpu::mgm::memfree_device(data.buffer_sptrs2);
        gpu::mgm::memfree_device(data.buffer_spmm);
    }
    d_applyg_x_cluster.clear();
    d_applyg_y_cluster.clear();
    d_applyg_xs_pointers.clear();
    d_applyg_ys_pointers.clear();
    d_applyg_n_dofs_interfaces.clear();
    d_applyg_D2Cs_pointers.clear();
    tm_gpumem.stop();

    tm_gpuhostmem.start();
    for(size_t d = 0; d < n_domains; d++)
    {
        per_domain_stuff & data = domain_data[d];
        data.h_L_sp.clear();
        data.h_LH_sp.clear();
        data.h_U_sp.clear();
        data.h_UH_sp.clear();
        data.h_Bperm_sp.clear();
        data.h_applyc_x.clear();
        data.h_applyc_y.clear();
    }
    tm_gpuhostmem.stop();

    tm_wait.start();
    gpu::mgm::device_wait();
    tm_wait.stop();

    tm_total.stop();

    stage = 0;

    print_timer("Destroy total", tm_total);
    print_timer("Destroy   descriptors", tm_descriptors);
    print_timer("Destroy   gpulibs", tm_gpulibs);
    print_timer("Destroy   gpumem", tm_gpumem);
    print_timer("Destroy     gpumempool", tm_gpumempool);
    print_timer("Destroy   gpuhostmem", tm_gpuhostmem);
    print_timer("Destroy   wait", tm_wait);
}

template <typename T, typename I>
void TotalFETIExplicitAcc<T,I>::info()
{
    // DualOperatorInfo sum, min, max;
    size_t minF = INT32_MAX, maxF = 0, sumF = 0;
    for (size_t d = 0; d < n_domains; ++d) {
        minF = std::min(minF, domain_data[d].d_F.nrows * domain_data[d].d_F.ncols * sizeof(T));
        maxF = std::max(maxF, domain_data[d].d_F.nrows * domain_data[d].d_F.ncols * sizeof(T));
        sumF += domain_data[d].d_F.nrows * domain_data[d].d_F.ncols * sizeof(T);
    }

//    TotalFETIImplicit<T>::reduceInfo(sum, min, max);
    Communication::allReduce(&minF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_MIN);
    Communication::allReduce(&maxF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_MAX);
    Communication::allReduce(&sumF, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_SUM);

    eslog::info(" = EXPLICIT TOTAL FETI OPERATOR ON GPU                                                       = \n");
//    TotalFETIImplicit<T>::printInfo(sum, min, max);
    config->ecfdescription->forEachParameters([](ECFParameter * param){
        std::string name = param->name;
        for(char & c : name) c = std::toupper(c);
        eslog::info(" =   %-50s       %+30s = \n", name.c_str(), param->getValue().c_str());
    });
    eslog::info(" =   F MEMORY [MB]                                            %8.2f <%8.2f - %8.2f> = \n", (double)sumF / n_domains / 1024. / 1024., minF / 1024. / 1024., maxF / 1024. / 1024.);
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}





template <typename T, typename I>
void TotalFETIExplicitAcc<T,I>::set(const step::Step &step)
{
    if(stage != 1) eslog::error("set: invalid order of operations in dualop\n");

    my_timer tm_total, tm_gpuinit, tm_gpuset, tm_vecresize, tm_sizecalc, tm_gpucreate, tm_mainloop_outer, tm_mainloop_inner, tm_Kreg_combine, tm_solver_commit, tm_fact_symbolic, tm_descriptors, tm_buffersize, tm_alloc, tm_alloc_host, tm_alloc_device, tm_setpointers, tm_Bperm, tm_get_factors, tm_extract, tm_transpose, tm_copyin, tm_kernels_preprocess, tm_trs1, tm_trs2, tm_gemm, tm_applystuff, tm_poolalloc, tm_wait;

    tm_total.start();
    n_domains = feti.K.size();

    if(config->queue_count == QUEUE_COUNT::PER_DOMAIN) n_queues = n_domains;
    if(config->queue_count == QUEUE_COUNT::PER_THREAD) n_queues = omp_get_max_threads();

    tm_gpuinit.start();
    gpu::mgm::init_gpu(device);
    tm_gpuinit.stop();

    tm_gpuset.start();
    gpu::mgm::set_device(device);
    tm_gpuset.stop();

    tm_vecresize.start();
    queues.resize(n_queues);
    handles_dense.resize(n_queues);
    handles_sparse.resize(n_queues);
    domain_data.resize(n_domains);
    d_Fs_allocated.resize(is_f_triangles_shared ? (n_domains - 1) / 2 + 1 : n_domains);
    tm_vecresize.stop();
    
    tm_sizecalc.start();
    for(size_t d = 0; d < n_domains; d++)
    {
        per_domain_stuff & data = domain_data[d];
        data.n_dofs_domain = feti.B1[d].ncols;
        data.n_dofs_interface = feti.B1[d].nrows;
        data.ld_domain = ((data.n_dofs_domain - 1) / align_elem + 1) * align_elem;
        data.ld_interface = ((data.n_dofs_interface - 1) / align_elem + 1) * align_elem;
        data.ld_X = (order_X == 'R' ? data.ld_interface : data.ld_domain);
    }
    if(is_f_triangles_shared)
    {
        std::vector<size_t> domain_idxs_sorted_by_f_size_desc(n_domains);
        for(size_t d = 0; d < n_domains; d++) domain_idxs_sorted_by_f_size_desc[d] = d;
        std::sort(domain_idxs_sorted_by_f_size_desc.rbegin(), domain_idxs_sorted_by_f_size_desc.rend(), [&](size_t dl, size_t dr){ return domain_data[dl].n_dofs_interface < domain_data[dr].n_dofs_interface; });
        for(size_t i = 0; i < n_domains; i++)
        {
            size_t d = domain_idxs_sorted_by_f_size_desc[i];
            size_t d_bigger = domain_idxs_sorted_by_f_size_desc[(i / 2) * 2];
            domain_data[d].allocated_F_index = i / 2;
            domain_data[d].hermitian_F_fill = (i % 2 == 0 ? 'U' : 'L');
            domain_data[d].ld_F = domain_data[d_bigger].ld_interface;
            domain_data[d].should_allocate_d_F = (i % 2 == 0);
        }
    }
    else
    {
        for(size_t d = 0; d < n_domains; d++)
        {
            domain_data[d].allocated_F_index = d;
            domain_data[d].hermitian_F_fill = 'U';
            domain_data[d].ld_F = domain_data[d].ld_interface;
            domain_data[d].should_allocate_d_F = true;
        }
    }
    {
        // approximate check if it is possible for the matrices to fit into memory
        size_t mem_needed = 0;
        for(size_t d = 0; d < n_domains; d++)
        {
            size_t mem_needed_f = Matrix_Dense<T,I>::memoryRequirement(domain_data[d].n_dofs_interface, domain_data[d].n_dofs_interface, domain_data[d].ld_F);
            if(is_f_triangles_shared) mem_needed_f /= 2;
            mem_needed += mem_needed_f;
        }
        size_t mem_capacity = gpu::mgm::get_device_memory_capacity();
        if(mem_needed > mem_capacity)
        {
            eslog::error("Not enough device memory for this input. Need %zu GiB, but have only %zu GiB\n", mem_needed >> 30, mem_capacity >> 30);
        }
    }
    tm_sizecalc.stop();

    // create gpu stream and libraries
    tm_gpucreate.start();
    gpu::mgm::queue_create(main_q);
    for(gpu::mgm::queue & q : queues) gpu::mgm::queue_create(q);
    for(size_t i = 0; i < n_queues; i++) gpu::dnblas::handle_create(handles_dense[i], queues[i]);
    for(size_t i = 0; i < n_queues; i++) gpu::spblas::handle_create(handles_sparse[i], queues[i]);
    tm_gpucreate.stop();

    tm_mainloop_outer.start();
    #pragma omp parallel for schedule(static,1) if(config->concurrency_set == CONCURRENCY::PARALLEL)
    for(size_t d = 0; d < n_domains; d++)
    {
        tm_mainloop_inner.start();

        I n_nz_factor;

        gpu::mgm::queue & q = queues[d % n_queues];
        gpu::dnblas::handle & hd = handles_dense[d % n_queues];
        gpu::spblas::handle & hs = handles_sparse[d % n_queues];
        per_domain_stuff & data = domain_data[d];

        gpu::spblas::descr_matrix_dense & descr_X = data.descr_X;
        gpu::spblas::descr_matrix_dense & descr_W = (is_factor1_sparse ? data.descr_Y : data.descr_X);
        gpu::spblas::descr_matrix_dense & descr_Z = (is_factor1_sparse == is_factor2_sparse ? data.descr_X : data.descr_Y);

        // Kreg = K + RegMat symbolic pattern
        tm_Kreg_combine.start();
        {
            math::combine(data.Kreg, feti.K[d], feti.RegMat[d]);
            if constexpr(utils::is_real<T>())    data.Kreg.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
            if constexpr(utils::is_complex<T>()) data.Kreg.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE;
            data.Kreg.shape = feti.K[d].shape;
        }
        tm_Kreg_combine.stop();

        // commit Kreg to solver (just symbolic pattern present now)
        tm_solver_commit.start();
        {
            data.solver_Kreg.commit(data.Kreg);
        }
        tm_solver_commit.stop();

        // symbolic factorization
        tm_fact_symbolic.start();
        {
            data.solver_Kreg.symbolicFactorization();
            n_nz_factor = data.solver_Kreg.getFactorNnz();
        }
        tm_fact_symbolic.stop();

        // create descriptors
        tm_descriptors.start();
        {
            if(trsm1_use_L)  gpu::spblas::descr_matrix_csr_create<T,I>(data.descr_L_sp,  data.n_dofs_domain, data.n_dofs_domain, n_nz_factor, 'L');
            if(trsm1_use_LH) gpu::spblas::descr_matrix_csr_create<T,I>(data.descr_LH_sp, data.n_dofs_domain, data.n_dofs_domain, n_nz_factor, 'U');
            if(trsm2_use_U)  gpu::spblas::descr_matrix_csr_create<T,I>(data.descr_U_sp,  data.n_dofs_domain, data.n_dofs_domain, n_nz_factor, 'U');
            if(trsm2_use_UH) gpu::spblas::descr_matrix_csr_create<T,I>(data.descr_UH_sp, data.n_dofs_domain, data.n_dofs_domain, n_nz_factor, 'L');
            gpu::spblas::descr_matrix_csr_create<T,I>(data.descr_Bperm_sp, data.n_dofs_interface, data.n_dofs_domain, feti.B1[d].nnz, 'N');
            if(trsm1_use_L  && is_factor1_dense) gpu::spblas::descr_matrix_dense_create<T,I>(data.descr_L_dn,  data.n_dofs_domain, data.n_dofs_domain, data.ld_domain, 'R');
            if(trsm1_use_LH && is_factor1_dense) gpu::spblas::descr_matrix_dense_create<T,I>(data.descr_LH_dn, data.n_dofs_domain, data.n_dofs_domain, data.ld_domain, 'R');
            if(trsm2_use_U  && is_factor2_dense) gpu::spblas::descr_matrix_dense_create<T,I>(data.descr_U_dn,  data.n_dofs_domain, data.n_dofs_domain, data.ld_domain, 'R');
            if(trsm2_use_UH && is_factor2_dense) gpu::spblas::descr_matrix_dense_create<T,I>(data.descr_UH_dn, data.n_dofs_domain, data.n_dofs_domain, data.ld_domain, 'R');
            if(true)   gpu::spblas::descr_matrix_dense_create<T,I>(data.descr_X, data.n_dofs_domain, data.n_dofs_interface, data.ld_X, order_X);
            if(need_Y) gpu::spblas::descr_matrix_dense_create<T,I>(data.descr_Y, data.n_dofs_domain, data.n_dofs_interface, data.ld_X, order_X);
            gpu::spblas::descr_matrix_dense_create<T,I>(data.descr_F, data.n_dofs_interface, data.n_dofs_interface, data.ld_F, order_F);
            if(is_factor1_sparse)                 gpu::spblas::descr_sparse_trsm_create(data.descr_sparse_trsm1);
            if(is_factor2_sparse && is_path_trsm) gpu::spblas::descr_sparse_trsm_create(data.descr_sparse_trsm2);
        }
        tm_descriptors.stop();

        // buffersize
        tm_buffersize.start();
        {
            std::vector<size_t> buffer_requirements;

            if(trsm1_use_L  && is_factor1_dense)                          gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_L_sp,  data.descr_L_dn,  buffer_requirements.emplace_back(), nullptr, 'B');
            if(trsm1_use_LH && is_factor1_dense && !can_use_LH_is_U_d_dn) gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_LH_sp, data.descr_LH_dn, buffer_requirements.emplace_back(), nullptr, 'B');
            if(trsm2_use_U  && is_factor2_dense)                          gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_U_sp,  data.descr_U_dn,  buffer_requirements.emplace_back(), nullptr, 'B');
            if(trsm2_use_UH && is_factor2_dense && !can_use_UH_is_L_d_dn) gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_UH_sp, data.descr_UH_dn, buffer_requirements.emplace_back(), nullptr, 'B');
            gpu::spblas::sparse_to_dense<T,I>(hs, 'T', data.descr_Bperm_sp, data.descr_X, buffer_requirements.emplace_back(), nullptr, 'B');

            if(trsm1_use_L  && is_factor1_sparse) gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_L_sp,  descr_X, descr_W, data.descr_sparse_trsm1, data.buffersize_sptrs1, data.buffer_sptrs1, 'B');
            if(trsm1_use_LH && is_factor1_sparse) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_LH_sp, descr_X, descr_W, data.descr_sparse_trsm1, data.buffersize_sptrs1, data.buffer_sptrs1, 'B');
            if(trsm2_use_U  && is_factor2_sparse) gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_U_sp,  descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersize_sptrs2, data.buffer_sptrs2, 'B');
            if(trsm2_use_UH && is_factor2_sparse) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_UH_sp, descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersize_sptrs2, data.buffer_sptrs2, 'B');
            if(is_path_trsm) gpu::spblas::mm<T,I>(hs, 'N', 'N', data.descr_Bperm_sp, descr_Z, data.descr_F, data.buffersize_spmm, data.buffer_spmm, 'B');

            gpu::dnblas::buffer_collect_size(hd, buffer_requirements.emplace_back(), [&](){
                T * dummyptrT = reinterpret_cast<T*>(sizeof(T));
                if(trsm1_use_L  && is_factor1_dense) gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, dummyptrT, data.ld_domain, 'R', 'N', 'L', dummyptrT, data.ld_X, order_X, 'N');
                if(trsm1_use_LH && is_factor1_dense) gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, dummyptrT, data.ld_domain, 'R', 'H', 'U', dummyptrT, data.ld_X, order_X, 'N');
                if(trsm2_use_U  && is_factor2_dense) gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, dummyptrT, data.ld_domain, 'R', 'N', 'U', dummyptrT, data.ld_X, order_X, 'N');
                if(trsm2_use_UH && is_factor2_dense) gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, dummyptrT, data.ld_domain, 'R', 'H', 'L', dummyptrT, data.ld_X, order_X, 'N');
                if(is_path_herk) gpu::dnblas::herk<T,I>(hd, data.n_dofs_interface, data.n_dofs_domain, dummyptrT, data.ld_X, order_X, 'H', dummyptrT, data.ld_F, order_F, data.hermitian_F_fill);
                if( is_system_hermitian) gpu::dnblas::hemv(hd, data.n_dofs_interface, dummyptrT, data.ld_F, order_F, 'N', data.hermitian_F_fill, dummyptrT, dummyptrT);
                if(!is_system_hermitian) gpu::dnblas::gemv(hd, data.n_dofs_interface, data.n_dofs_interface, dummyptrT, data.ld_F, order_F, 'N', dummyptrT, dummyptrT);
            });

            gpu::mgm::queue_wait(q);

            data.buffersize_other = *std::max_element(buffer_requirements.begin(), buffer_requirements.end());
        }
        tm_buffersize.stop();

        // permanent allocations for the lifetime of the object
        tm_alloc.start();
        {
            // host pinned memory
            tm_alloc_host.start();
            if(trsm1_use_L || solver_get_L)           data.h_L_sp .resize(data.n_dofs_domain, data.n_dofs_domain, n_nz_factor);
            if(trsm2_use_U || solver_get_U)           data.h_U_sp .resize(data.n_dofs_domain, data.n_dofs_domain, n_nz_factor);
            if(trsm1_use_LH && !can_use_LH_is_U_h_sp) data.h_LH_sp.resize(data.n_dofs_domain, data.n_dofs_domain, n_nz_factor);
            if(trsm2_use_UH && !can_use_UH_is_L_h_sp) data.h_UH_sp.resize(data.n_dofs_domain, data.n_dofs_domain, n_nz_factor);
            if(                 can_use_LH_is_U_h_sp) data.h_LH_sp.shallowCopy(data.h_U_sp);
            if(                 can_use_UH_is_L_h_sp) data.h_UH_sp.shallowCopy(data.h_L_sp);
            data.h_Bperm_sp.resize(data.n_dofs_interface, data.n_dofs_domain, feti.B1[d].nnz);
            if(config->apply_scatter_gather_where == DEVICE::CPU) data.h_applyc_x.resize(data.n_dofs_interface);
            if(config->apply_scatter_gather_where == DEVICE::CPU) data.h_applyc_y.resize(data.n_dofs_interface);
            tm_alloc_host.stop();

            // device memory
            tm_alloc_device.start();
            if(trsm1_use_L)                           data.d_L_sp .resize(data.n_dofs_domain, data.n_dofs_domain, n_nz_factor);
            if(trsm2_use_U)                           data.d_U_sp .resize(data.n_dofs_domain, data.n_dofs_domain, n_nz_factor);
            if(trsm1_use_LH && !can_use_LH_is_U_d_sp) data.d_LH_sp.resize(data.n_dofs_domain, data.n_dofs_domain, n_nz_factor);
            if(trsm2_use_UH && !can_use_UH_is_L_d_sp) data.d_UH_sp.resize(data.n_dofs_domain, data.n_dofs_domain, n_nz_factor);
            if(trsm1_use_LH &&  can_use_LH_is_U_d_sp) data.d_LH_sp.shallowCopy(data.d_U_sp);
            if(trsm2_use_UH &&  can_use_UH_is_L_d_sp) data.d_UH_sp.shallowCopy(data.d_L_sp);
            data.d_Bperm_sp.resize(data.n_dofs_interface, data.n_dofs_domain, feti.B1[d].nnz);
            data.d_apply_x.resize(data.n_dofs_interface);
            data.d_apply_y.resize(data.n_dofs_interface);
            data.d_applyg_D2C.resize(data.n_dofs_interface);
            if(is_factor1_sparse)                 data.buffer_sptrs1 = gpu::mgm::memalloc_device(data.buffersize_sptrs1);
            if(is_factor2_sparse && is_path_trsm) data.buffer_sptrs2 = gpu::mgm::memalloc_device(data.buffersize_sptrs2);
            if(                     is_path_trsm) data.buffer_spmm   = gpu::mgm::memalloc_device(data.buffersize_spmm);
            if(data.should_allocate_d_F) d_Fs_allocated[data.allocated_F_index].resize(data.n_dofs_interface + 1, data.n_dofs_interface, data.ld_F);
            tm_alloc_device.stop();
        }
        tm_alloc.stop();

        // set the pointers inside the descriptors of some matrices
        tm_setpointers.start();
        {
            if(trsm1_use_L)  gpu::spblas::descr_matrix_csr_link_data(data.descr_L_sp,  data.d_L_sp);
            if(trsm1_use_LH) gpu::spblas::descr_matrix_csr_link_data(data.descr_LH_sp, data.d_LH_sp);
            if(trsm2_use_U)  gpu::spblas::descr_matrix_csr_link_data(data.descr_U_sp,  data.d_U_sp);
            if(trsm2_use_UH) gpu::spblas::descr_matrix_csr_link_data(data.descr_UH_sp, data.d_UH_sp);
            gpu::spblas::descr_matrix_csr_link_data(data.descr_Bperm_sp, data.d_Bperm_sp);
        }
        tm_setpointers.stop();

        // prepare matrices on host
        tm_Bperm.start();
        {
            Permutation<I> perm;
            perm.resize(data.n_dofs_domain);
            data.solver_Kreg.getPermutation(perm);
            math::permuteColumns(data.h_Bperm_sp, feti.B1[d], perm);
        }
        tm_Bperm.stop();

        // extract symbolic pattern from the factor
        tm_get_factors.start();
        {
            tm_extract.start();
            if(solver_get_L) data.solver_Kreg.getFactorL(data.h_L_sp, true, false);
            if(solver_get_U) data.solver_Kreg.getFactorU(data.h_U_sp, true, false);
            tm_extract.stop();
            tm_transpose.start();
            if(need_conjtrans_L2LH)
            {
                data.transmap_L2LH.resize(data.h_L_sp.nnz);
                math::conjTransposeMapSetup(data.h_LH_sp, data.transmap_L2LH, data.h_L_sp);
            }
            if(need_conjtrans_U2UH)
            {
                data.transmap_U2UH.resize(data.h_U_sp.nnz);
                math::conjTransposeMapSetup(data.h_UH_sp, data.transmap_U2UH, data.h_U_sp);
            }

            tm_transpose.stop();
        }
        tm_get_factors.stop();

        // copy some matrices to device
        tm_copyin.start();
        {
            if(trsm1_use_L)                           gpu::mgm::copy_submit_h2d(q, data.d_L_sp,  data.h_L_sp,  true, false);
            if(trsm2_use_U)                           gpu::mgm::copy_submit_h2d(q, data.d_U_sp,  data.h_U_sp,  true, false);
            if(trsm1_use_LH && !can_use_LH_is_U_d_sp) gpu::mgm::copy_submit_h2d(q, data.d_LH_sp, data.h_LH_sp, true, false);
            if(trsm2_use_UH && !can_use_UH_is_L_d_sp) gpu::mgm::copy_submit_h2d(q, data.d_UH_sp, data.h_UH_sp, true, false);
            gpu::mgm::copy_submit_h2d(q, data.d_Bperm_sp, data.h_Bperm_sp, true, true);
            if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        }
        tm_copyin.stop();
        
        // proprocessing stage of the kernels
        tm_kernels_preprocess.start();
        {
            tm_trs1.start();
            if(trsm1_use_L  && is_factor1_sparse) gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_L_sp,  descr_X, descr_W, data.descr_sparse_trsm1, data.buffersize_sptrs1, data.buffer_sptrs1, 'P');
            if(trsm1_use_LH && is_factor1_sparse) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_LH_sp, descr_X, descr_W, data.descr_sparse_trsm1, data.buffersize_sptrs1, data.buffer_sptrs1, 'P');
            if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
            tm_trs1.stop();

            if(is_path_trsm)
            {
                tm_trs2.start();
                if(trsm2_use_U  && is_factor2_sparse) gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_U_sp,  descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersize_sptrs2, data.buffer_sptrs2, 'P');
                if(trsm2_use_UH && is_factor2_sparse) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_UH_sp, descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersize_sptrs2, data.buffer_sptrs2, 'P');
                if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                tm_trs2.stop();
                
                tm_gemm.start();
                gpu::spblas::mm<T,I>(hs, 'N', 'N', data.descr_Bperm_sp, descr_Z, data.descr_F, data.buffersize_spmm, data.buffer_spmm, 'P');
                if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                tm_gemm.stop();
            }
        }
        tm_kernels_preprocess.stop();

        tm_mainloop_inner.stop();
    }
    tm_mainloop_outer.stop();

    for(size_t d = 0; d < n_domains; d++)
    {
        per_domain_stuff & data = domain_data[d];
        data.d_F.shallowCopy(d_Fs_allocated[data.allocated_F_index]);
        data.d_F.nrows = data.n_dofs_interface;
        data.d_F.ncols = data.n_dofs_interface;
        if((data.hermitian_F_fill == 'L') == (order_F == 'R')) data.d_F.vals += data.d_F.get_ld();
    }

    // clean up the mess from buggy openmp in clang
    utils::run_dummy_parallel_region();
    
    // some stuff needed for apply
    tm_applystuff.start();
    if(config->apply_scatter_gather_where == DEVICE::GPU) d_applyg_x_cluster.resize(feti.lambdas.size);
    if(config->apply_scatter_gather_where == DEVICE::GPU) d_applyg_y_cluster.resize(feti.lambdas.size);
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
            h_applyg_xs_pointers.vals[d] = domain_data[d].d_apply_x.vals;
            h_applyg_ys_pointers.vals[d] = domain_data[d].d_apply_y.vals;
            h_applyg_n_dofs_interfaces.vals[d] = domain_data[d].n_dofs_interface;
            h_applyg_D2Cs_pointers.vals[d] = domain_data[d].d_applyg_D2C.vals;
        }
        gpu::mgm::copy_submit_h2d(main_q, d_applyg_xs_pointers,       h_applyg_xs_pointers);
        gpu::mgm::copy_submit_h2d(main_q, d_applyg_ys_pointers,       h_applyg_ys_pointers);
        gpu::mgm::copy_submit_h2d(main_q, d_applyg_n_dofs_interfaces, h_applyg_n_dofs_interfaces);
        gpu::mgm::copy_submit_h2d(main_q, d_applyg_D2Cs_pointers,     h_applyg_D2Cs_pointers);
        for(size_t d = 0; d < n_domains; d++) gpu::mgm::copy_submit_h2d(main_q, domain_data[d].d_applyg_D2C.vals, feti.D2C[d].data(), feti.D2C[d].size());
        if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
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
                per_domain_stuff & data = domain_data[d];
                size_t mem_pool_size_request_domain = 0;
                if(is_factor1_dense && trsm1_use_L)                           mem_pool_size_request_domain += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
                if(is_factor2_dense && trsm2_use_U)                           mem_pool_size_request_domain += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
                if(is_factor1_dense && trsm1_use_LH && !can_use_LH_is_U_d_dn) mem_pool_size_request_domain += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
                if(is_factor2_dense && trsm2_use_UH && !can_use_UH_is_L_d_dn) mem_pool_size_request_domain += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
                if(true)   mem_pool_size_request_domain += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_interface, data.ld_X);
                if(need_Y) mem_pool_size_request_domain += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_interface, data.ld_X);
                if(need_f_tmp) Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_interface, data.n_dofs_interface, data.ld_F);
                mem_pool_size_request_domain += domain_data[d].buffersize_other;
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
    print_timer("Set       sizecalc", tm_sizecalc);
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

    my_timer tm_total, tm_mainloop_outer, tm_mainloop_inner, tm_Kreg_combine, tm_solver_commit, tm_fact_numeric, tm_get_factors, tm_extract, tm_transpose, tm_allocinpool, tm_setpointers, tm_copyin, tm_descr_update, tm_sp2dn, tm_kernels_compute, tm_trs1, tm_trs2, tm_gemm, tm_fcopy, tm_syrk, tm_freeinpool, tm_freeinpool_exec, tm_compute_d, tm_wait;

    tm_total.start();

    gpu::mgm::set_device(device);

    tm_mainloop_outer.start();
    #pragma omp parallel for schedule(static,1) if(config->concurrency_update == CONCURRENCY::PARALLEL)
    for(size_t d = 0; d < n_domains; d++)
    {
        tm_mainloop_inner.start();

        gpu::mgm::queue & q = queues[d % n_queues];
        gpu::dnblas::handle & hd = handles_dense[d % n_queues];
        gpu::spblas::handle & hs = handles_sparse[d % n_queues];
        per_domain_stuff & data = domain_data[d];

        std::unique_ptr<Matrix_Dense<T,I,cbmba_d>> & d_X = data.d_X;
        std::unique_ptr<Matrix_Dense<T,I,cbmba_d>> & d_W = (is_factor1_sparse ? data.d_Y : data.d_X);
        // std::unique_ptr<Matrix_Dense<T,I,cbmba_d>> & d_Z = (is_factor1_sparse == is_factor2_sparse ? data.d_X : data.d_Y);
        gpu::spblas::descr_matrix_dense & descr_X = data.descr_X;
        gpu::spblas::descr_matrix_dense & descr_W = (is_factor1_sparse ? data.descr_Y : data.descr_X);
        gpu::spblas::descr_matrix_dense & descr_Z = (is_factor1_sparse == is_factor2_sparse ? data.descr_X : data.descr_Y);

        void * buffer_other = nullptr;

        // Kreg = K + RegMat numeric values
        tm_Kreg_combine.start();
        {
            math::sumCombined(data.Kreg, T{1.0}, feti.K[d], feti.RegMat[d]);
        }
        tm_Kreg_combine.stop();

        // commit Kreg to solver (with numeric values)
        tm_solver_commit.start();
        {
            data.solver_Kreg.commit(data.Kreg);
        }
        tm_solver_commit.stop();
        
        // numeric factorization
        tm_fact_numeric.start();
        {
            data.solver_Kreg.numericalFactorization();
        }
        tm_fact_numeric.stop();

        // extract values from numeric factor
        tm_get_factors.start();
        {
            tm_extract.start();
            if(solver_get_L) data.solver_Kreg.getFactorL(data.h_L_sp, false, true);
            if(solver_get_U) data.solver_Kreg.getFactorU(data.h_U_sp, false, true);
            tm_extract.stop();
            tm_transpose.start();
            if(need_conjtrans_L2LH) math::conjTransposeMapUse(data.h_LH_sp, data.transmap_L2LH, data.h_L_sp);
            if(need_conjtrans_U2UH) math::conjTransposeMapUse(data.h_UH_sp, data.transmap_U2UH, data.h_U_sp);
            tm_transpose.stop();
        }
        tm_get_factors.stop();

        // temporary allocations using the memory pool
        tm_allocinpool.start();
        {
            cbmba_d ator_d(*cbmba_res_device, align_B);
            cbmba_res_device->do_transaction([&](){
                if(is_factor1_dense && trsm1_use_L)  data.d_L_dn  = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
                if(is_factor2_dense && trsm2_use_U)  data.d_U_dn  = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
                if(is_factor1_dense && trsm1_use_LH) data.d_LH_dn = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
                if(is_factor2_dense && trsm2_use_UH) data.d_UH_dn = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
                if(true)   data.d_X = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
                if(need_Y) data.d_Y = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
                data.d_F_tmp = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);

                if(is_factor1_dense && trsm1_use_L)                           data.d_L_dn ->resize(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
                if(is_factor2_dense && trsm2_use_U)                           data.d_U_dn ->resize(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
                if(is_factor1_dense && trsm1_use_LH && !can_use_LH_is_U_d_dn) data.d_LH_dn->resize(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
                if(is_factor2_dense && trsm2_use_UH && !can_use_UH_is_L_d_dn) data.d_UH_dn->resize(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
                if(is_factor1_dense && trsm1_use_LH &&  can_use_LH_is_U_d_dn) data.d_LH_dn->shallowCopy(*data.d_U_dn);
                if(is_factor2_dense && trsm2_use_UH &&  can_use_UH_is_L_d_dn) data.d_UH_dn->shallowCopy(*data.d_L_dn);
                if(order_X == 'R')           data.d_X->resize(data.n_dofs_domain,    data.n_dofs_interface, data.ld_X);
                if(order_X == 'C')           data.d_X->resize(data.n_dofs_interface, data.n_dofs_domain,    data.ld_X);
                if(order_X == 'R' && need_Y) data.d_Y->resize(data.n_dofs_domain,    data.n_dofs_interface, data.ld_X);
                if(order_X == 'C' && need_Y) data.d_Y->resize(data.n_dofs_interface, data.n_dofs_domain,    data.ld_X);
                if( need_f_tmp) data.d_F_tmp->resize(data.n_dofs_interface, data.n_dofs_interface, data.ld_F);
                if(!need_f_tmp) data.d_F_tmp->shallowCopy(data.d_F);
                buffer_other = cbmba_res_device->allocate(data.buffersize_other, align_B);
            });
        }
        tm_allocinpool.stop();

        gpu::dnblas::buffer_set(hd, buffer_other, data.buffersize_other);

        // set the pointers inside the descriptors of the rest of the matrices
        tm_setpointers.start();
        {
            if(trsm1_use_L  && is_factor1_dense) gpu::spblas::descr_matrix_dense_link_data(data.descr_L_dn,  *data.d_L_dn);
            if(trsm1_use_LH && is_factor1_dense) gpu::spblas::descr_matrix_dense_link_data(data.descr_LH_dn, *data.d_LH_dn);
            if(trsm2_use_U  && is_factor2_dense) gpu::spblas::descr_matrix_dense_link_data(data.descr_U_dn,  *data.d_U_dn);
            if(trsm2_use_UH && is_factor2_dense) gpu::spblas::descr_matrix_dense_link_data(data.descr_UH_dn, *data.d_UH_dn);
            if(true)   gpu::spblas::descr_matrix_dense_link_data(data.descr_X, *data.d_X);
            if(need_Y) gpu::spblas::descr_matrix_dense_link_data(data.descr_Y, *data.d_Y);
            if(is_path_trsm) gpu::spblas::descr_matrix_dense_link_data(data.descr_F, *data.d_F_tmp);
        }
        tm_setpointers.stop();

        // copy the new factors to device
        tm_copyin.start();
        {
            if(trsm1_use_L)                           gpu::mgm::copy_submit_h2d(q, data.d_L_sp,  data.h_L_sp,  false, true);
            if(trsm2_use_U)                           gpu::mgm::copy_submit_h2d(q, data.d_U_sp,  data.h_U_sp,  false, true);
            if(trsm1_use_LH && !can_use_LH_is_U_d_sp) gpu::mgm::copy_submit_h2d(q, data.d_LH_sp, data.h_LH_sp, false, true);
            if(trsm2_use_UH && !can_use_UH_is_L_d_sp) gpu::mgm::copy_submit_h2d(q, data.d_UH_sp, data.h_UH_sp, false, true);
            if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        }
        tm_copyin.stop();

        // update sparse trsm descriptors to reflect the new matrix values, possibly re-preprocess
        tm_descr_update.start();
        {
            if(trsm1_use_L  && is_factor1_sparse) gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_L_sp,  descr_X, descr_W, data.descr_sparse_trsm1, data.buffersize_sptrs1, data.buffer_sptrs1, 'U');
            if(trsm1_use_LH && is_factor1_sparse) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_LH_sp, descr_X, descr_W, data.descr_sparse_trsm1, data.buffersize_sptrs1, data.buffer_sptrs1, 'U');
            if(trsm2_use_U  && is_factor2_sparse) gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_U_sp,  descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersize_sptrs2, data.buffer_sptrs2, 'U');
            if(trsm2_use_UH && is_factor2_sparse) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_UH_sp, descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersize_sptrs2, data.buffer_sptrs2, 'U');
            if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        }
        tm_descr_update.stop();

        // sparse to dense on device
        tm_sp2dn.start();
        {
            if(trsm1_use_L  && is_factor1_dense)                          gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_L_sp,  data.descr_L_dn,  data.buffersize_other, buffer_other, 'C');
            if(trsm2_use_U  && is_factor2_dense)                          gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_U_sp,  data.descr_U_dn,  data.buffersize_other, buffer_other, 'C');
            if(trsm1_use_LH && is_factor1_dense && !can_use_LH_is_U_d_dn) gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_LH_sp, data.descr_LH_dn, data.buffersize_other, buffer_other, 'C');
            if(trsm2_use_UH && is_factor2_dense && !can_use_UH_is_L_d_dn) gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_UH_sp, data.descr_UH_dn, data.buffersize_other, buffer_other, 'C');
            gpu::spblas::sparse_to_dense<T,I>(hs, 'T', data.descr_Bperm_sp, data.descr_X, data.buffersize_other, buffer_other, 'C');
            if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        }
        tm_sp2dn.stop();

        // perform the actual assembly
        tm_kernels_compute.start();
        {
            tm_trs1.start();
            if(is_factor1_sparse && trsm1_use_L)  gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_L_sp,  descr_X, descr_W, data.descr_sparse_trsm1, data.buffersize_sptrs1, data.buffer_sptrs1, 'C');
            if(is_factor1_sparse && trsm1_use_LH) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_LH_sp, descr_X, descr_W, data.descr_sparse_trsm1, data.buffersize_sptrs1, data.buffer_sptrs1, 'C');
            if(is_factor1_dense  && trsm1_use_L)  gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, data.d_L_dn->vals,  data.ld_domain, 'R', 'N', 'L', d_X->vals, data.ld_X, order_X, 'N');
            if(is_factor1_dense  && trsm1_use_LH) gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, data.d_LH_dn->vals, data.ld_domain, 'R', 'H', 'U', d_X->vals, data.ld_X, order_X, 'N');
            if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
            tm_trs1.stop();

            if(is_path_herk)
            {
                tm_syrk.start();
                gpu::dnblas::herk<T,I>(hd, data.n_dofs_interface, data.n_dofs_domain, d_W->vals, data.ld_X, order_X, 'H', data.d_F_tmp->vals, data.ld_F, order_F, data.hermitian_F_fill);
                if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                tm_syrk.stop();
            }
            if(is_path_trsm)
            {
                tm_trs2.start();
                if(is_factor2_sparse && trsm2_use_U)  gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_U_sp,  descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersize_sptrs2, data.buffer_sptrs2, 'C');
                if(is_factor2_sparse && trsm2_use_UH) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_UH_sp, descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersize_sptrs2, data.buffer_sptrs2, 'C');
                if(is_factor2_dense  && trsm2_use_U)  gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, data.d_U_dn->vals,  data.ld_domain, 'R', 'N', 'U', d_W->vals, data.ld_X, order_X, 'N');
                if(is_factor2_dense  && trsm2_use_UH) gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, data.d_UH_dn->vals, data.ld_domain, 'R', 'H', 'L', d_W->vals, data.ld_X, order_X, 'N');
                if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                tm_trs2.stop();

                tm_gemm.start();
                gpu::spblas::mm<T,I>(hs, 'N', 'N', data.descr_Bperm_sp, descr_Z, data.descr_F, data.buffersize_spmm, data.buffer_spmm, 'C');
                if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                tm_gemm.stop();

                tm_fcopy.start();
                if(need_f_tmp) gpu::kernels::copy_matrix_triangle(q, data.d_F, *data.d_F_tmp, data.hermitian_F_fill, order_F);
                if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                tm_fcopy.stop();
            }
        }
        tm_kernels_compute.stop();

        gpu::dnblas::buffer_unset(hd);

        // free the temporary memory from the pool
        tm_freeinpool.start();
        gpu::mgm::submit_host_function(q, [&,d,buffer_other](){
            tm_freeinpool_exec.start();
            data.d_L_dn.reset();
            data.d_LH_dn.reset();
            data.d_U_dn.reset();
            data.d_UH_dn.reset();
            data.d_X.reset();
            data.d_Y.reset();
            data.d_F_tmp.reset();
            cbmba_res_device->deallocate(buffer_other);
            tm_freeinpool_exec.stop();
        });
        if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
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
        #pragma omp parallel for schedule(static,1) if(config->concurrency_update == CONCURRENCY::PARALLEL)
        for(size_t d = 0; d < n_domains; d++) {
            domain_data[d].solver_Kreg.solve(feti.f[d], Kplus_fs[d]);
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
    print_timer("Update          fcopy", tm_fcopy);
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

    if(config->apply_scatter_gather_where == DEVICE::GPU)
    {
        my_timer tm_total, tm_copyin, tm_scatter, tm_mv_outer, tm_mv, tm_zerofill, tm_gather, tm_copyout, tm_wait;

        tm_total.start();

        // copy x_cluster to device
        tm_copyin.start();
        gpu::mgm::copy_submit_h2d(main_q, d_applyg_x_cluster, x_cluster);
        if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
        tm_copyin.stop();

        // scatter
        tm_scatter.start();
        gpu::kernels::DCmap_scatter(main_q, d_applyg_xs_pointers, d_applyg_n_dofs_interfaces, d_applyg_x_cluster, d_applyg_D2Cs_pointers);
        if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
        tm_scatter.stop();

        // apply
        tm_mv_outer.start();
        gpu::mgm::queue_async_barrier({main_q}, queues);
        #pragma omp parallel for schedule(static,1) if(config->concurrency_apply == CONCURRENCY::PARALLEL)
        for(size_t d = 0; d < n_domains; d++)
        {
            gpu::mgm::queue & q = queues[d % n_queues];
            gpu::dnblas::handle & hd = handles_dense[d % n_queues];
            per_domain_stuff & data = domain_data[d];

            tm_mv.start();
            if( is_system_hermitian) gpu::dnblas::hemv<T,I>(hd, data.d_F.nrows, data.d_F.vals, data.d_F.get_ld(), order_F, 'N', data.hermitian_F_fill, data.d_apply_x.vals, data.d_apply_y.vals);
            if(!is_system_hermitian) gpu::dnblas::gemv<T,I>(hd, data.d_F.nrows, data.d_F.ncols, data.d_F.vals, data.d_F.get_ld(), order_F, 'N', data.d_apply_x.vals, data.d_apply_y.vals);
            if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
            tm_mv.stop();
        }
        tm_mv_outer.stop();

        // zerofill y_cluster on device
        tm_zerofill.start();
        gpu::mgm::memset_submit(main_q, d_applyg_y_cluster.vals, d_applyg_y_cluster.size * sizeof(T), 0);
        if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
        tm_zerofill.stop();

        gpu::mgm::queue_async_barrier(queues, {main_q});

        // gather
        tm_gather.start();
        gpu::kernels::DCmap_gather(main_q, d_applyg_ys_pointers, d_applyg_n_dofs_interfaces, d_applyg_y_cluster, d_applyg_D2Cs_pointers);
        if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
        tm_gather.stop();

        // copy y_cluster from device
        tm_copyout.start();
        gpu::mgm::copy_submit_d2h(main_q, y_cluster, d_applyg_y_cluster);
        if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
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
    
    if(config->apply_scatter_gather_where == DEVICE::CPU)
    {
        my_timer tm_total, tm_apply_outer, tm_apply_inner, tm_scatter, tm_copyin, tm_mv, tm_copyout, tm_zerofill, tm_wait, tm_gather, tm_gather_inner;

        tm_total.start();

        tm_apply_outer.start();
        #pragma omp parallel for schedule(static,1) if(config->concurrency_apply == CONCURRENCY::PARALLEL)
        for(size_t d = 0; d < n_domains; d++)
        {
            tm_apply_inner.start();

            gpu::mgm::queue & q = queues[d % n_queues];
            gpu::dnblas::handle & hd = handles_dense[d % n_queues];
            per_domain_stuff & data = domain_data[d];

            tm_scatter.start();
            for(I i = 0; i < data.n_dofs_interface; i++) data.h_applyc_x.vals[i] = x_cluster.vals[feti.D2C[d][i]];
            tm_scatter.stop();

            tm_copyin.start();
            gpu::mgm::copy_submit_h2d(q, data.d_apply_x, data.h_applyc_x);
            if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
            tm_copyin.stop();

            tm_mv.start();
            if( is_system_hermitian) gpu::dnblas::hemv<T,I>(hd, data.d_F.nrows, data.d_F.vals, data.d_F.get_ld(), order_F, 'N', data.hermitian_F_fill, data.d_apply_x.vals, data.d_apply_y.vals);
            if(!is_system_hermitian) gpu::dnblas::gemv<T,I>(hd, data.d_F.nrows, data.d_F.ncols, data.d_F.vals, data.d_F.get_ld(), order_F, 'N', data.d_apply_x.vals, data.d_apply_y.vals);
            if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
            tm_mv.stop();

            tm_copyout.start();
            gpu::mgm::copy_submit_d2h(q, data.h_applyc_y, data.d_apply_y);
            if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
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
        #pragma omp parallel for schedule(static,1) if(config->concurrency_apply == CONCURRENCY::PARALLEL)
        for(size_t d = 0; d < n_domains; d++)
        {
            tm_gather_inner.start();
            per_domain_stuff & data = domain_data[d];
            for(I i = 0; i < data.n_dofs_interface; i++)
            {
                if constexpr(utils::is_real<T>())
                {
                    #pragma omp atomic
                    y_cluster.vals[feti.D2C[d][i]] += data.h_applyc_y.vals[i];
                }
                if constexpr(utils::is_complex<T>())
                {
                    #pragma omp atomic
                    utils::real_ref(y_cluster.vals[feti.D2C[d][i]]) += utils::real_ref(data.h_applyc_y.vals[i]);
                    #pragma omp atomic
                    utils::imag_ref(y_cluster.vals[feti.D2C[d][i]]) += utils::imag_ref(data.h_applyc_y.vals[i]);
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
    #pragma omp parallel for schedule(static,1) if(config->concurrency_set == CONCURRENCY::PARALLEL)
    for (size_t d = 0; d < feti.K.size(); ++d) {
        Vector_Dense<T,I> z;
        z.resize(y[d]);
        applyBt(feti, d, x, z, T{-1});
        math::add(z, T{1}, feti.f[d]);
        domain_data[d].solver_Kreg.solve(z, y[d]);
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



template <typename T, typename I>
template<typename U>
void TotalFETIExplicitAcc<T,I>::replace_if_default(U & val, U replace_with)
{
    if(val == U::DEFAULT) val = replace_with;
}



template <typename T, typename I>
void TotalFETIExplicitAcc<T,I>::config_replace_defaults()
{
    if(gpu::mgm::get_implementation() == gpu::mgm::gpu_wrapper_impl::CUDA)
    {
        replace_if_default(config->trsm1_factor_storage,       MATRIX_STORAGE::SPARSE);
        replace_if_default(config->trsm2_factor_storage,       MATRIX_STORAGE::SPARSE);
        replace_if_default(config->trsm1_solve_type,           TRSM1_SOLVE_TYPE::LHH);
        replace_if_default(config->trsm2_solve_type,           TRSM2_SOLVE_TYPE::U);
        replace_if_default(config->trsm_rhs_sol_order,         MATRIX_ORDER::ROW_MAJOR);
        replace_if_default(config->path_if_hermitian,          PATH_IF_HERMITIAN::HERK);
        replace_if_default(config->f_sharing_if_hermitian,     TRIANGLE_MATRIX_SHARING::SHARED);
        replace_if_default(config->queue_count,                QUEUE_COUNT::PER_THREAD);
        replace_if_default(config->apply_scatter_gather_where, DEVICE::GPU);
        
        bool is_spsm_used = (config->trsm1_factor_storage == MATRIX_STORAGE::SPARSE || (config->trsm2_factor_storage == MATRIX_STORAGE::SPARSE && config->path_if_hermitian == PATH_IF_HERMITIAN::TRSM));
        bool is_spsm_concurrent_working;
        if(gpu::spblas::get_implementation() == gpu::spblas::spblas_wrapper_impl::CUSPARSE_LEGACY) {
            is_spsm_concurrent_working = true;
        }
        else if(gpu::spblas::get_implementation() == gpu::spblas::spblas_wrapper_impl::CUSPARSE_MODERN) {
            is_spsm_concurrent_working = false;
        }
        else {
            eslog::error("Unexpected gpu sparse blas implementation\n");
        }
        bool need_seq = (is_spsm_used && !is_spsm_concurrent_working);
        CONCURRENCY concurrency_set_update = (need_seq ? CONCURRENCY::SEQ_WAIT : CONCURRENCY::PARALLEL);
        replace_if_default(config->concurrency_set,            concurrency_set_update);
        replace_if_default(config->concurrency_update,         concurrency_set_update);
        replace_if_default(config->concurrency_apply,          CONCURRENCY::SEQ_CONTINUE);
    }
    if(gpu::mgm::get_implementation() == gpu::mgm::gpu_wrapper_impl::ROCM)
    {
        replace_if_default(config->concurrency_set,            CONCURRENCY::PARALLEL);
        replace_if_default(config->concurrency_update,         CONCURRENCY::PARALLEL);
        replace_if_default(config->concurrency_apply,          CONCURRENCY::SEQ_CONTINUE);
        replace_if_default(config->trsm1_factor_storage,       MATRIX_STORAGE::SPARSE);
        replace_if_default(config->trsm2_factor_storage,       MATRIX_STORAGE::SPARSE);
        replace_if_default(config->trsm1_solve_type,           TRSM1_SOLVE_TYPE::LHH);
        replace_if_default(config->trsm2_solve_type,           TRSM2_SOLVE_TYPE::U);
        replace_if_default(config->trsm_rhs_sol_order,         MATRIX_ORDER::ROW_MAJOR);
        replace_if_default(config->path_if_hermitian,          PATH_IF_HERMITIAN::HERK);
        replace_if_default(config->f_sharing_if_hermitian,     TRIANGLE_MATRIX_SHARING::SHARED);
        replace_if_default(config->queue_count,                QUEUE_COUNT::PER_THREAD);
        replace_if_default(config->apply_scatter_gather_where, DEVICE::GPU);
    }
}



#define INSTANTIATE_T_I(T,I) \
template class TotalFETIExplicitAcc<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

        // INSTANTIATE_T(float)
        INSTANTIATE_T(double)
        // INSTANTIATE_T(std::complex<float>)
        // INSTANTIATE_T(std::complex<double>)

    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I

}
