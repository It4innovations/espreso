
#include <numeric>
#include <limits>

#include "totalfeti.gpu.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "basis/utilities/minmaxavg.h"

#include "my_timer.h"

namespace espreso {

using CONCURRENCY             = DualOperatorGpuConfig::CONCURRENCY;
using MATRIX_STORAGE          = DualOperatorGpuConfig::MATRIX_STORAGE;
using TRS1_SOLVE_TYPE         = DualOperatorGpuConfig::TRS1_SOLVE_TYPE;
using TRS2_SOLVE_TYPE         = DualOperatorGpuConfig::TRS2_SOLVE_TYPE;
using MATRIX_ORDER            = DualOperatorGpuConfig::MATRIX_ORDER;
using PATH_IF_HERMITIAN       = DualOperatorGpuConfig::PATH_IF_HERMITIAN;
using TRIANGLE_MATRIX_SHARING = DualOperatorGpuConfig::TRIANGLE_MATRIX_SHARING;
using QUEUE_COUNT             = DualOperatorGpuConfig::QUEUE_COUNT;
using DEVICE                  = DualOperatorGpuConfig::DEVICE;
using TIMERS                  = DualOperatorGpuConfig::TIMERS;
using MEMORY_INFO             = DualOperatorGpuConfig::MEMORY_INFO;

template <typename T, typename I>
TotalFETIGpu<T,I>::TotalFETIGpu(FETI<T> &feti, DualOperatorStrategy strategy)
: DualOperator<T>(feti), n_domains(0), n_queues(0), mem_pool_device(nullptr)
{
    if(stage != 0) eslog::error("init: invalid order of operations in dualop\n");
    
    is_explicit = (strategy == DualOperatorStrategy::EXPLICIT);
    is_implicit = (strategy == DualOperatorStrategy::IMPLICIT);

    config = &feti.configuration.dual_operator_gpu_config;
    config_replace_defaults();
    if(config->trsm_rhs_sol_order == MATRIX_ORDER::ROW_MAJOR) order_X = 'R';
    if(config->trsm_rhs_sol_order == MATRIX_ORDER::COL_MAJOR) order_X = 'C';
    order_F = order_X;
    is_system_hermitian = (DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_LOWER || DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_UPPER);
    is_factor1_dense = (config->trs1_factor_storage == MATRIX_STORAGE::DENSE);
    is_factor2_dense = (config->trs2_factor_storage == MATRIX_STORAGE::DENSE);
    is_factor1_sparse = (config->trs1_factor_storage == MATRIX_STORAGE::SPARSE);
    is_factor2_sparse = (config->trs2_factor_storage == MATRIX_STORAGE::SPARSE);
    is_path_trsm = (config->path_if_hermitian == PATH_IF_HERMITIAN::TRSM);
    is_path_herk = (config->path_if_hermitian == PATH_IF_HERMITIAN::HERK);
    do_herk = (is_path_herk && is_explicit);

    bool is_sptrsm_inplace = (gpu::spblas::get_place_trsm() == gpu::spblas::place::IN_PLACE);
    // bool is_sptrsm_outofplace = (gpu::spblas::get_place_trsm() == gpu::spblas::place::OUT_OF_PLACE);
    is_trsm1_inplace = (is_factor1_dense || (is_factor1_sparse && is_sptrsm_inplace));
    is_trsm2_inplace = (is_factor2_dense || (is_factor2_sparse && is_sptrsm_inplace));
    is_trsm1_outofplace = !is_trsm1_inplace; // or, equivalently, is_factor1_sparse && is_sptrsm_outofplace
    is_trsm2_outofplace = !is_trsm2_inplace; // or, equivalently, is_factor2_sparse && is_sptrsm_outofplace

    bool trs1_use_L  = (config->trs1_solve_type == TRS1_SOLVE_TYPE::L);
    bool trs1_use_LH = (config->trs1_solve_type == TRS1_SOLVE_TYPE::LHH);
    bool trs2_use_U  = (config->trs2_solve_type == TRS2_SOLVE_TYPE::U   && (is_path_trsm || is_implicit));
    bool trs2_use_UH = (config->trs2_solve_type == TRS2_SOLVE_TYPE::UHH && (is_path_trsm || is_implicit));
    need_X = is_explicit;
    need_Y = (is_explicit && (is_trsm1_outofplace || is_trsm2_outofplace));
    need_F = is_explicit;
    Solver_Factors sym = DirectSparseSolver<T,I>::factorsSymmetry();
    solver_get_L = (sym == Solver_Factors::HERMITIAN_LOWER || sym == Solver_Factors::NONSYMMETRIC_BOTH);
    solver_get_U = (sym == Solver_Factors::HERMITIAN_UPPER || sym == Solver_Factors::NONSYMMETRIC_BOTH);
    bool conjtrans_on_host   = (config->transpose_where == DEVICE::CPU);
    bool conjtrans_on_device = (config->transpose_where == DEVICE::GPU);
    is_f_triangles_shared = (is_system_hermitian && config->f_sharing_if_hermitian == TRIANGLE_MATRIX_SHARING::SHARED);
    need_f_tmp = (need_F && is_f_triangles_shared && is_path_trsm);
    timers_basic = (config->timers == TIMERS::BASIC || config->timers == TIMERS::ALL);
    timers_detailed = (config->timers == TIMERS::ALL);
    memory_info_basic = (config->memory_info == MEMORY_INFO::BASIC || config->memory_info == MEMORY_INFO::ALL);
    memory_info_detailed = (config->memory_info == MEMORY_INFO::ALL);

    bool need_conjtrans_L2LH = (is_system_hermitian && solver_get_L && (trs1_use_LH || trs2_use_U)) || (!is_system_hermitian && trs1_use_LH);
    bool need_conjtrans_U2UH = (is_system_hermitian && solver_get_U && (trs2_use_UH || trs1_use_L)) || (!is_system_hermitian && trs2_use_UH);
    do_conjtrans_L2LH_h = (need_conjtrans_L2LH && conjtrans_on_host);
    do_conjtrans_U2UH_h = (need_conjtrans_U2UH && conjtrans_on_host);
    do_conjtrans_L2LH_d = (need_conjtrans_L2LH && conjtrans_on_device);
    do_conjtrans_U2UH_d = (need_conjtrans_U2UH && conjtrans_on_device);

    do_trsm1_sp = (is_factor1_sparse && is_explicit);
    do_trsm2_sp = (is_factor2_sparse && is_explicit && is_path_trsm);

    do_trsv1_sp = (is_factor1_sparse && is_implicit);
    do_trsv2_sp = (is_factor2_sparse && is_implicit);

    bool do_trs1_sp_L  = (trs1_use_L  && is_factor1_sparse);
    bool do_trs1_sp_LH = (trs1_use_LH && is_factor1_sparse);
    bool do_trs2_sp_U  = (trs2_use_U  && is_factor2_sparse);
    bool do_trs2_sp_UH = (trs2_use_UH && is_factor2_sparse);
    bool do_trs1_dn_L  = (trs1_use_L  && is_factor1_dense);
    bool do_trs1_dn_LH = (trs1_use_LH && is_factor1_dense);
    bool do_trs2_dn_U  = (trs2_use_U  && is_factor2_dense);
    bool do_trs2_dn_UH = (trs2_use_UH && is_factor2_dense);

    do_trsm1_sp_L  = (do_trs1_sp_L  && is_explicit);
    do_trsm1_sp_LH = (do_trs1_sp_LH && is_explicit);
    do_trsm2_sp_U  = (do_trs2_sp_U  && is_explicit);
    do_trsm2_sp_UH = (do_trs2_sp_UH && is_explicit);
    do_trsm1_dn_L  = (do_trs1_dn_L  && is_explicit);
    do_trsm1_dn_LH = (do_trs1_dn_LH && is_explicit);
    do_trsm2_dn_U  = (do_trs2_dn_U  && is_explicit);
    do_trsm2_dn_UH = (do_trs2_dn_UH && is_explicit);

    do_trsv1_sp_L  = (do_trs1_sp_L  && is_implicit);
    do_trsv1_sp_LH = (do_trs1_sp_LH && is_implicit);
    do_trsv2_sp_U  = (do_trs2_sp_U  && is_implicit);
    do_trsv2_sp_UH = (do_trs2_sp_UH && is_implicit);
    do_trsv1_dn_L  = (do_trs1_dn_L  && is_implicit);
    do_trsv1_dn_LH = (do_trs1_dn_LH && is_implicit);
    do_trsv2_dn_U  = (do_trs2_dn_U  && is_implicit);
    do_trsv2_dn_UH = (do_trs2_dn_UH && is_implicit);

    do_trs1_sp = (is_factor1_sparse);
    do_trs2_sp = (is_factor2_sparse && (is_path_trsm || is_implicit));
    do_mm      = (is_path_trsm  && is_explicit);

    do_apply_hemv = (is_explicit &&  is_system_hermitian);
    do_apply_gemv = (is_explicit && !is_system_hermitian);

    do_alloc_h_sp_L   = (solver_get_L || (trs1_use_L && conjtrans_on_host));
    do_alloc_h_sp_U   = (solver_get_U || (trs2_use_U && conjtrans_on_host));
    bool can_use_LH_is_U_h_sp = (is_system_hermitian && do_alloc_h_sp_U);
    bool can_use_UH_is_L_h_sp = (is_system_hermitian && do_alloc_h_sp_L);
    do_alloc_h_sp_LH  = (!can_use_LH_is_U_h_sp && trs1_use_LH && conjtrans_on_host);
    do_alloc_h_sp_UH  = (!can_use_UH_is_L_h_sp && trs2_use_UH && conjtrans_on_host);
    do_link_h_sp_LH_U = ( can_use_LH_is_U_h_sp);
    do_link_h_sp_UH_L = ( can_use_UH_is_L_h_sp);
    bool is_present_h_sp_L  = (do_alloc_h_sp_L);
    bool is_present_h_sp_U  = (do_alloc_h_sp_U);
    bool is_present_h_sp_LH = (do_alloc_h_sp_LH || do_link_h_sp_LH_U);
    bool is_present_h_sp_UH = (do_alloc_h_sp_UH || do_link_h_sp_UH_L);

    do_alloc_d_sp_L = (trs1_use_L || ((trs1_use_LH || (trs2_use_U && is_system_hermitian)) && solver_get_L && conjtrans_on_device));
    do_alloc_d_sp_U = (trs2_use_U || ((trs2_use_UH || (trs1_use_L && is_system_hermitian)) && solver_get_U && conjtrans_on_device));
    bool can_use_LH_is_U_d_sp = (is_system_hermitian && do_alloc_d_sp_U);
    bool can_use_UH_is_L_d_sp = (is_system_hermitian && do_alloc_d_sp_L);
    do_alloc_d_sp_LH  = (!can_use_LH_is_U_d_sp && trs1_use_LH);
    do_alloc_d_sp_UH  = (!can_use_UH_is_L_d_sp && trs2_use_UH);
    do_link_d_sp_LH_U = ( can_use_LH_is_U_d_sp);
    do_link_d_sp_UH_L = ( can_use_UH_is_L_d_sp);
    bool is_present_d_sp_L  = (do_alloc_d_sp_L);
    bool is_present_d_sp_U  = (do_alloc_d_sp_U);
    bool is_present_d_sp_LH = (do_alloc_d_sp_LH || do_link_d_sp_LH_U);
    bool is_present_d_sp_UH = (do_alloc_d_sp_UH || do_link_d_sp_UH_L);

    do_alloc_d_dn_L  = (is_factor1_dense && trs1_use_L);
    do_alloc_d_dn_U  = (is_factor2_dense && trs2_use_U);
    bool can_use_LH_is_U_d_dn = (is_system_hermitian && do_alloc_d_dn_U);
    bool can_use_UH_is_L_d_dn = (is_system_hermitian && do_alloc_d_dn_L);
    do_alloc_d_dn_LH  = (!can_use_LH_is_U_d_dn && is_factor1_dense && trs1_use_LH);
    do_alloc_d_dn_UH  = (!can_use_UH_is_L_d_dn && is_factor2_dense && trs2_use_UH);
    do_link_d_dn_LH_U = ( can_use_LH_is_U_d_dn);
    do_link_d_dn_UH_L = ( can_use_UH_is_L_d_dn);

    do_sp2dn_L  = (trs1_use_L  && is_factor1_dense);
    do_sp2dn_U  = (trs2_use_U  && is_factor2_dense);
    do_sp2dn_LH = (trs1_use_LH && is_factor1_dense && !(can_use_LH_is_U_d_dn && do_sp2dn_U));
    do_sp2dn_UH = (trs2_use_UH && is_factor2_dense && !(can_use_UH_is_L_d_dn && do_sp2dn_L));
    do_sp2dn_X  = need_X;

    do_descr_sp_L  = (do_sp2dn_L  || do_trs1_sp_L  || do_conjtrans_L2LH_d);
    do_descr_sp_LH = (do_sp2dn_LH || do_trs1_sp_LH || do_conjtrans_L2LH_d);
    do_descr_sp_U  = (do_sp2dn_U  || do_trs2_sp_U  || do_conjtrans_U2UH_d);
    do_descr_sp_UH = (do_sp2dn_UH || do_trs2_sp_UH || do_conjtrans_U2UH_d);
    do_descr_dn_L  = (do_sp2dn_L);
    do_descr_dn_LH = (do_sp2dn_LH);
    do_descr_dn_U  = (do_sp2dn_U);
    do_descr_dn_UH = (do_sp2dn_UH);

    do_copyin_L = (is_present_h_sp_L && is_present_d_sp_L);
    do_copyin_U = (is_present_h_sp_U && is_present_d_sp_U);
    do_copyin_LH = (is_present_h_sp_LH && is_present_d_sp_LH && !(can_use_LH_is_U_d_sp && do_copyin_U));
    do_copyin_UH = (is_present_h_sp_UH && is_present_d_sp_UH && !(can_use_UH_is_L_d_sp && do_copyin_L));

    if(is_path_herk && !is_system_hermitian) eslog::error("cannot do herk path with non-hermitian system\n");

    device = gpu::mgm::get_device_by_mpi(info::mpi::rank, info::mpi::size);

    stage = 1;
}

template <typename T, typename I>
TotalFETIGpu<T,I>::~TotalFETIGpu()
{
    gpu::mgm::set_device(device);

    my_timer tm_total(timers_basic);
    my_timer tm_descriptors(timers_detailed), tm_gpulibs(timers_detailed), tm_gpumem(timers_detailed), tm_gpumempool(timers_detailed), tm_gpuhostmem(timers_detailed), tm_wait(timers_detailed);

    tm_total.start();

    tm_descriptors.start();
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        gpu::spblas::handle & hs = handles_sparse[di % n_queues];

        gpu::spblas::descr_matrix_csr_destroy(hs, data.descr_L_sp);
        gpu::spblas::descr_matrix_csr_destroy(hs, data.descr_LH_sp);
        gpu::spblas::descr_matrix_csr_destroy(hs, data.descr_U_sp);
        gpu::spblas::descr_matrix_csr_destroy(hs, data.descr_UH_sp);
        gpu::spblas::descr_matrix_csr_destroy(hs, data.descr_Bperm_sp);
        gpu::spblas::descr_matrix_dense_destroy(hs, data.descr_L_dn);
        gpu::spblas::descr_matrix_dense_destroy(hs, data.descr_LH_dn);
        gpu::spblas::descr_matrix_dense_destroy(hs, data.descr_U_dn);
        gpu::spblas::descr_matrix_dense_destroy(hs, data.descr_UH_dn);
        gpu::spblas::descr_matrix_dense_destroy(hs, data.descr_X);
        gpu::spblas::descr_matrix_dense_destroy(hs, data.descr_Y);
        gpu::spblas::descr_matrix_dense_destroy(hs, data.descr_F);
        gpu::spblas::descr_matrix_dense_destroy(hs, data.descr_F_tmp);
        gpu::spblas::descr_sparse_trsm_destroy(hs, data.descr_sparse_trsm1);
        gpu::spblas::descr_sparse_trsm_destroy(hs, data.descr_sparse_trsm2);
        gpu::spblas::descr_vector_dense_destroy(hs, data.descr_xvec);
        gpu::spblas::descr_vector_dense_destroy(hs, data.descr_yvec);
        gpu::spblas::descr_vector_dense_destroy(hs, data.descr_zvec);
        gpu::spblas::descr_vector_dense_destroy(hs, data.descr_wvec);
        gpu::spblas::descr_sparse_trsv_destroy(hs, data.descr_sparse_trsv1);
        gpu::spblas::descr_sparse_trsv_destroy(hs, data.descr_sparse_trsv2);
        gpu::spblas::descr_sparse_mv_destroy(hs, data.descr_sparse_spmv1);
        gpu::spblas::descr_sparse_mv_destroy(hs, data.descr_sparse_spmv2);
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
    for(auto & f : d_Fs_allocated) {
        f.clear();
    }
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        data.d_L_sp.clear();
        data.d_LH_sp.clear();
        data.d_U_sp.clear();
        data.d_UH_sp.clear();
        data.d_Bperm_sp.clear();
        data.d_applyg_D2C.clear();
        data.d_apply_x.clear();
        data.d_apply_y.clear();
        gpu::mgm::memfree_device(data.buffer_sptrsm1_persistent);
        gpu::mgm::memfree_device(data.buffer_sptrsm2_persistent);
        gpu::mgm::memfree_device(data.buffer_sptrsv1_persistent);
        gpu::mgm::memfree_device(data.buffer_sptrsv2_persistent);
        gpu::mgm::memfree_device(data.buffer_spmm);
        gpu::mgm::memfree_device(data.buffer_transL2LH);
        gpu::mgm::memfree_device(data.buffer_transU2UH);
        gpu::mgm::memfree_device(data.buffer_spmv1);
        gpu::mgm::memfree_device(data.buffer_spmv2);
    }
    d_applyg_x_cluster.clear();
    d_applyg_y_cluster.clear();
    d_applyg_xs_pointers.clear();
    d_applyg_ys_pointers.clear();
    d_applyg_n_dofs_interfaces.clear();
    d_applyg_D2Cs_pointers.clear();
    tm_gpumem.stop();

    tm_gpuhostmem.start();
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
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
void TotalFETIGpu<T,I>::info()
{
    if(stage < 2) eslog::error("info: invalid order of operations in dualop\n");

    eslog::info(" = EXPLICIT TOTAL FETI OPERATOR ON GPU                                                       = \n");
    config->ecfdescription->forEachParameters([](ECFParameter * param){
        std::string name = param->name;
        for(char & c : name) c = std::toupper(c);
        eslog::info(" =   %-50s       %+30s = \n", name.c_str(), param->getValue().c_str());
    });
    eslog::info(minmaxavg<double>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.d_F.nrows * data.d_F.get_ld() * sizeof(T) / (1024.0 * 1024.0); }).to_string("  F MEMORY [MB]").c_str());
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_domain; }).to_string("  Domain volume [dofs]").c_str());
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_interface; }).to_string("  Domain surface [dofs]").c_str());
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_nz_factor; }).to_string("  Factor nnz").c_str());
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}





template <typename T, typename I>
void TotalFETIGpu<T,I>::set(const step::Step &step)
{
    if(stage != 1) eslog::error("set: invalid order of operations in dualop\n");

    my_timer tm_total(timers_basic), tm_gpuinit(timers_basic), tm_gpuset(timers_basic), tm_vecresize(timers_basic), tm_sizecalc(timers_basic), tm_gpucreate(timers_basic), tm_libinit(timers_basic), tm_mainloop_outer(timers_basic), tm_applystuff(timers_basic), tm_mempool_reqs(timers_basic), tm_poolalloc(timers_basic), tm_preprocess(timers_basic), tm_wait(timers_basic);
    my_timer tm_mainloop_inner(timers_detailed), tm_Kreg_combine(timers_detailed), tm_solver_commit(timers_detailed), tm_fact_symbolic(timers_detailed), tm_descriptors(timers_detailed), tm_buffersize(timers_detailed), tm_alloc(timers_detailed), tm_alloc_host(timers_detailed), tm_alloc_device(timers_detailed), tm_setpointers(timers_detailed), tm_Bperm(timers_detailed), tm_get_factors(timers_detailed), tm_extract(timers_detailed), tm_trans_cpu(timers_detailed), tm_copyin(timers_detailed), tm_prepr_poolalloc(timers_detailed), tm_prepr_work(timers_detailed), tm_prepr_allocinpool(timers_detailed), tm_prepr_kernels(timers_detailed), tm_trans_gpu(timers_detailed), tm_trsm1(timers_detailed), tm_trsm2(timers_detailed), tm_gemm(timers_detailed), tm_spmv1(timers_detailed), tm_spmv2(timers_detailed), tm_trsv1(timers_detailed), tm_trsv2(timers_detailed), tm_prepr_pooldelete(timers_detailed);

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
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_stuff & data = domain_data[di];
        data.n_dofs_domain = feti.B1[di].ncols;
        data.n_dofs_interface = feti.B1[di].nrows;
        data.ld_domain = ((data.n_dofs_domain - 1) / align_elem + 1) * align_elem;
        data.ld_interface = ((data.n_dofs_interface - 1) / align_elem + 1) * align_elem;
        data.ld_X = (order_X == 'R' ? data.ld_interface : data.ld_domain);
    }
    if(need_F) {
        if(is_f_triangles_shared) {
            std::vector<size_t> domain_idxs_sorted_by_f_size_desc(n_domains);
            for(size_t di = 0; di < n_domains; di++) {
                domain_idxs_sorted_by_f_size_desc[di] = di;
            }
            std::sort(domain_idxs_sorted_by_f_size_desc.rbegin(), domain_idxs_sorted_by_f_size_desc.rend(), [&](size_t dl, size_t dr){ return domain_data[dl].n_dofs_interface < domain_data[dr].n_dofs_interface; });
            for(size_t i = 0; i < n_domains; i++) {
                size_t di = domain_idxs_sorted_by_f_size_desc[i];
                size_t di_bigger = domain_idxs_sorted_by_f_size_desc[(i / 2) * 2];
                domain_data[di].allocated_F_index = i / 2;
                domain_data[di].hermitian_F_fill = (i % 2 == 0 ? 'U' : 'L');
                domain_data[di].ld_F = domain_data[di_bigger].ld_interface;
                domain_data[di].should_allocate_d_F = (i % 2 == 0);
            }
        }
        else {
            for(size_t di = 0; di < n_domains; di++) {
                domain_data[di].allocated_F_index = di;
                domain_data[di].hermitian_F_fill = 'U';
                domain_data[di].ld_F = domain_data[di].ld_interface;
                domain_data[di].should_allocate_d_F = true;
            }
        }
        {
            // approximate check if it is possible for the matrices to fit into memory
            size_t mem_needed = 0;
            for(size_t di = 0; di < n_domains; di++) {
                size_t mem_needed_f = Matrix_Dense<T,I>::memoryRequirement(domain_data[di].n_dofs_interface, domain_data[di].n_dofs_interface, domain_data[di].ld_F);
                if(is_f_triangles_shared) mem_needed_f /= 2;
                mem_needed += mem_needed_f;
            }
            size_t mem_capacity = gpu::mgm::get_device_memory_capacity();
            if(mem_needed > mem_capacity) {
                eslog::error("Not enough device memory for this input. Need %zu GiB, but have only %zu GiB\n", mem_needed >> 30, mem_capacity >> 30);
            }
        }
    }
    else
    {
        for(size_t di = 0; di < n_domains; di++) {
            domain_data[di].allocated_F_index = -1;
            domain_data[di].hermitian_F_fill = '_';
            domain_data[di].ld_F = 0;
            domain_data[di].should_allocate_d_F = false;
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

    tm_libinit.start();
    gpu::dnblas::init_library(main_q);
    gpu::mgm::queue_wait(main_q);
    tm_libinit.stop();

    tm_mainloop_outer.start();
    #pragma omp parallel for schedule(static,1) if(config->concurrency_set == CONCURRENCY::PARALLEL)
    for(size_t di = 0; di < n_domains; di++) {
        tm_mainloop_inner.start();

        gpu::mgm::queue & q = queues[di % n_queues];
        gpu::dnblas::handle & hd = handles_dense[di % n_queues];
        gpu::spblas::handle & hs = handles_sparse[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        gpu::spblas::descr_matrix_dense & descr_X = data.descr_X;
        gpu::spblas::descr_matrix_dense & descr_W = (is_trsm1_outofplace ? data.descr_Y : data.descr_X);
        gpu::spblas::descr_matrix_dense & descr_Z = (is_trsm1_outofplace == is_trsm2_outofplace ? data.descr_X : data.descr_Y);

        // Kreg = K + RegMat symbolic pattern
        tm_Kreg_combine.start();
        {
            math::combine(data.Kreg, feti.K[di], feti.RegMat[di]);
            if constexpr(utils::is_real<T>())    data.Kreg.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
            if constexpr(utils::is_complex<T>()) data.Kreg.type = Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE;
            data.Kreg.shape = feti.K[di].shape;
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
            data.n_nz_factor = data.solver_Kreg.getFactorNnz();
        }
        tm_fact_symbolic.stop();

        // create descriptors
        tm_descriptors.start();
        {
            gpu::spblas::handle & hs = handles_sparse[di % n_queues];

            if(do_descr_sp_L)  gpu::spblas::descr_matrix_csr_create<T,I>(hs, data.descr_L_sp,  data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor, 'L');
            if(do_descr_sp_LH) gpu::spblas::descr_matrix_csr_create<T,I>(hs, data.descr_LH_sp, data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor, 'U');
            if(do_descr_sp_U)  gpu::spblas::descr_matrix_csr_create<T,I>(hs, data.descr_U_sp,  data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor, 'U');
            if(do_descr_sp_UH) gpu::spblas::descr_matrix_csr_create<T,I>(hs, data.descr_UH_sp, data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor, 'L');
            gpu::spblas::descr_matrix_csr_create<T,I>(hs, data.descr_Bperm_sp, data.n_dofs_interface, data.n_dofs_domain, feti.B1[di].nnz, 'N');
            if(do_descr_dn_L)  gpu::spblas::descr_matrix_dense_create<T,I>(hs, data.descr_L_dn,  data.n_dofs_domain, data.n_dofs_domain, data.ld_domain, 'R');
            if(do_descr_dn_LH) gpu::spblas::descr_matrix_dense_create<T,I>(hs, data.descr_LH_dn, data.n_dofs_domain, data.n_dofs_domain, data.ld_domain, 'R');
            if(do_descr_dn_U)  gpu::spblas::descr_matrix_dense_create<T,I>(hs, data.descr_U_dn,  data.n_dofs_domain, data.n_dofs_domain, data.ld_domain, 'R');
            if(do_descr_dn_UH) gpu::spblas::descr_matrix_dense_create<T,I>(hs, data.descr_UH_dn, data.n_dofs_domain, data.n_dofs_domain, data.ld_domain, 'R');
            if(need_X) gpu::spblas::descr_matrix_dense_create<T,I>(hs, data.descr_X, data.n_dofs_domain, data.n_dofs_interface, data.ld_X, order_X);
            if(need_Y) gpu::spblas::descr_matrix_dense_create<T,I>(hs, data.descr_Y, data.n_dofs_domain, data.n_dofs_interface, data.ld_X, order_X);
            if(need_F)     gpu::spblas::descr_matrix_dense_create<T,I>(hs, data.descr_F,     data.n_dofs_interface, data.n_dofs_interface, data.ld_F, order_F);
            if(need_f_tmp) gpu::spblas::descr_matrix_dense_create<T,I>(hs, data.descr_F_tmp, data.n_dofs_interface, data.n_dofs_interface, data.ld_F, order_F);
            if(do_trsm1_sp) gpu::spblas::descr_sparse_trsm_create(hs, data.descr_sparse_trsm1);
            if(do_trsm2_sp) gpu::spblas::descr_sparse_trsm_create(hs, data.descr_sparse_trsm2);
            if(is_implicit) gpu::spblas::descr_vector_dense_create<T,I>(hs, data.descr_xvec, data.n_dofs_interface);
            if(is_implicit) gpu::spblas::descr_vector_dense_create<T,I>(hs, data.descr_yvec, data.n_dofs_interface);
            if(is_implicit) gpu::spblas::descr_vector_dense_create<T,I>(hs, data.descr_zvec, data.n_dofs_domain);
            if(is_implicit) gpu::spblas::descr_vector_dense_create<T,I>(hs, data.descr_wvec, data.n_dofs_domain);
            if(is_implicit) gpu::spblas::descr_sparse_trsv_create(hs, data.descr_sparse_trsv1);
            if(is_implicit) gpu::spblas::descr_sparse_trsv_create(hs, data.descr_sparse_trsv2);
            if(is_implicit) gpu::spblas::descr_sparse_mv_create(hs, data.descr_sparse_spmv1);
            if(is_implicit) gpu::spblas::descr_sparse_mv_create(hs, data.descr_sparse_spmv2);
        }
        tm_descriptors.stop();

        // buffersize
        tm_buffersize.start();
        {
            std::vector<size_t> buffersize_tmp_requirements_set;
            std::vector<size_t> buffersize_tmp_requirements_update;
            std::vector<size_t> buffersize_tmp_requirements_apply;
            buffersize_tmp_requirements_set.push_back(0);
            buffersize_tmp_requirements_update.push_back(0);
            buffersize_tmp_requirements_apply.push_back(0);

            if(do_conjtrans_L2LH_d) gpu::spblas::transpose<T,I>(hs, data.descr_LH_sp, data.descr_L_sp, true, data.buffersize_transL2LH, data.buffer_transL2LH, 'B');
            if(do_conjtrans_U2UH_d) gpu::spblas::transpose<T,I>(hs, data.descr_UH_sp, data.descr_U_sp, true, data.buffersize_transU2UH, data.buffer_transU2UH, 'B');

            {
                std::vector<size_t> * the_vector = nullptr;
                if(is_explicit) the_vector = &buffersize_tmp_requirements_update;
                if(is_implicit) the_vector = &buffersize_tmp_requirements_apply;
                if(do_sp2dn_L)  gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_L_sp,  data.descr_L_dn,  the_vector->emplace_back(), nullptr, 'B');
                if(do_sp2dn_U)  gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_U_sp,  data.descr_U_dn,  the_vector->emplace_back(), nullptr, 'B');
                if(do_sp2dn_LH) gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_LH_sp, data.descr_LH_dn, the_vector->emplace_back(), nullptr, 'B');
                if(do_sp2dn_UH) gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_UH_sp, data.descr_UH_dn, the_vector->emplace_back(), nullptr, 'B');
                if(do_sp2dn_X)  gpu::spblas::sparse_to_dense<T,I>(hs, 'T', data.descr_Bperm_sp, data.descr_X,  the_vector->emplace_back(), nullptr, 'B');
            }

            if(do_mm &&  need_f_tmp) gpu::spblas::mm<T,I>(hs, 'N', 'N', data.descr_Bperm_sp, descr_Z, data.descr_F_tmp, data.buffersize_spmm, data.buffer_spmm, 'B');
            if(do_mm && !need_f_tmp) gpu::spblas::mm<T,I>(hs, 'N', 'N', data.descr_Bperm_sp, descr_Z, data.descr_F,     data.buffersize_spmm, data.buffer_spmm, 'B');
            if(do_trsm1_sp_L)  gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_L_sp,  descr_X, descr_W, data.descr_sparse_trsm1, data.buffersizes_sptrsm1, nullptr, nullptr, 'B');
            if(do_trsm1_sp_LH) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_LH_sp, descr_X, descr_W, data.descr_sparse_trsm1, data.buffersizes_sptrsm1, nullptr, nullptr, 'B');
            if(do_trsm2_sp_U)  gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_U_sp,  descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersizes_sptrsm2, nullptr, nullptr, 'B');
            if(do_trsm2_sp_UH) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_UH_sp, descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersizes_sptrsm2, nullptr, nullptr, 'B');
            buffersize_tmp_requirements_set.push_back(data.buffersizes_sptrsm1.tmp_preprocess);
            buffersize_tmp_requirements_set.push_back(data.buffersizes_sptrsm2.tmp_preprocess);
            buffersize_tmp_requirements_update.push_back(data.buffersizes_sptrsm1.tmp_update);
            buffersize_tmp_requirements_update.push_back(data.buffersizes_sptrsm2.tmp_update);
            buffersize_tmp_requirements_update.push_back(data.buffersizes_sptrsm1.tmp_compute);
            buffersize_tmp_requirements_update.push_back(data.buffersizes_sptrsm2.tmp_compute);

            if(is_implicit) gpu::spblas::mv<T,I>(hs, 'T', data.descr_Bperm_sp, data.descr_xvec, data.descr_wvec, data.descr_sparse_spmv1, data.buffersize_spmv1, data.buffer_spmv1, 'B');
            if(is_implicit) gpu::spblas::mv<T,I>(hs, 'N', data.descr_Bperm_sp, data.descr_wvec, data.descr_yvec, data.descr_sparse_spmv2, data.buffersize_spmv2, data.buffer_spmv2, 'B');
            if(do_trsv1_sp_L)  gpu::spblas::trsv<T,I>(hs, 'N', data.descr_L_sp,  data.descr_wvec, data.descr_zvec, data.descr_sparse_trsv1, data.buffersizes_sptrsv1, nullptr, nullptr, 'B');
            if(do_trsv1_sp_LH) gpu::spblas::trsv<T,I>(hs, 'H', data.descr_LH_sp, data.descr_wvec, data.descr_zvec, data.descr_sparse_trsv1, data.buffersizes_sptrsv1, nullptr, nullptr, 'B');
            if(do_trsv2_sp_U)  gpu::spblas::trsv<T,I>(hs, 'N', data.descr_U_sp,  data.descr_zvec, data.descr_wvec, data.descr_sparse_trsv2, data.buffersizes_sptrsv2, nullptr, nullptr, 'B');
            if(do_trsv2_sp_UH) gpu::spblas::trsv<T,I>(hs, 'H', data.descr_UH_sp, data.descr_zvec, data.descr_wvec, data.descr_sparse_trsv2, data.buffersizes_sptrsv2, nullptr, nullptr, 'B');
            buffersize_tmp_requirements_set.push_back(data.buffersizes_sptrsv1.tmp_preprocess);
            buffersize_tmp_requirements_set.push_back(data.buffersizes_sptrsv2.tmp_preprocess);
            buffersize_tmp_requirements_update.push_back(data.buffersizes_sptrsv1.tmp_update);
            buffersize_tmp_requirements_update.push_back(data.buffersizes_sptrsv2.tmp_update);
            buffersize_tmp_requirements_apply.push_back(data.buffersizes_sptrsv1.tmp_compute);
            buffersize_tmp_requirements_apply.push_back(data.buffersizes_sptrsv2.tmp_compute);

            gpu::dnblas::buffer_collect_size(hd, buffersize_tmp_requirements_update.emplace_back(), [&](){
                T * dummyptrT = reinterpret_cast<T*>(sizeof(T));
                if(do_trsm1_dn_L)  gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, dummyptrT, data.ld_domain, 'R', 'N', 'L', dummyptrT, data.ld_X, order_X, 'N');
                if(do_trsm1_dn_LH) gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, dummyptrT, data.ld_domain, 'R', 'H', 'U', dummyptrT, data.ld_X, order_X, 'N');
                if(do_trsm2_dn_U)  gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, dummyptrT, data.ld_domain, 'R', 'N', 'U', dummyptrT, data.ld_X, order_X, 'N');
                if(do_trsm2_dn_UH) gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, dummyptrT, data.ld_domain, 'R', 'H', 'L', dummyptrT, data.ld_X, order_X, 'N');
                if(do_herk) gpu::dnblas::herk<T,I>(hd, data.n_dofs_interface, data.n_dofs_domain, dummyptrT, data.ld_X, order_X, 'H', dummyptrT, data.ld_F, order_F, data.hermitian_F_fill);
            });
            gpu::dnblas::buffer_collect_size(hd, buffersize_tmp_requirements_apply.emplace_back(), [&](){
                T * dummyptrT = reinterpret_cast<T*>(sizeof(T));
                if(do_apply_hemv) gpu::dnblas::hemv(hd, data.n_dofs_interface, dummyptrT, data.ld_F, order_F, 'N', data.hermitian_F_fill, dummyptrT, dummyptrT);
                if(do_apply_gemv) gpu::dnblas::gemv(hd, data.n_dofs_interface, data.n_dofs_interface, dummyptrT, data.ld_F, order_F, 'N', dummyptrT, dummyptrT);
                if(do_trsv1_dn_L)  gpu::dnblas::trsv<T,I>(hd, data.n_dofs_domain, dummyptrT, data.ld_domain, 'R', 'N', 'L', dummyptrT);
                if(do_trsv1_dn_LH) gpu::dnblas::trsv<T,I>(hd, data.n_dofs_domain, dummyptrT, data.ld_domain, 'R', 'H', 'U', dummyptrT);
                if(do_trsv2_dn_U)  gpu::dnblas::trsv<T,I>(hd, data.n_dofs_domain, dummyptrT, data.ld_domain, 'R', 'N', 'U', dummyptrT);
                if(do_trsv2_dn_UH) gpu::dnblas::trsv<T,I>(hd, data.n_dofs_domain, dummyptrT, data.ld_domain, 'R', 'H', 'L', dummyptrT);
            });

            gpu::mgm::queue_wait(q);

            data.buffersize_tmp_max_set    = *std::max_element(buffersize_tmp_requirements_set.begin(),    buffersize_tmp_requirements_set.end());
            data.buffersize_tmp_max_update = *std::max_element(buffersize_tmp_requirements_update.begin(), buffersize_tmp_requirements_update.end());
            data.buffersize_tmp_max_apply  = *std::max_element(buffersize_tmp_requirements_apply.begin(),  buffersize_tmp_requirements_apply.end());
        }
        tm_buffersize.stop();

        // permanent allocations for the lifetime of the object
        tm_alloc.start();
        {
            // host pinned memory
            tm_alloc_host.start();
            if(do_alloc_h_sp_L)   data.h_L_sp .resize(data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor);
            if(do_alloc_h_sp_U)   data.h_U_sp .resize(data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor);
            if(do_alloc_h_sp_LH)  data.h_LH_sp.resize(data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor);
            if(do_alloc_h_sp_UH)  data.h_UH_sp.resize(data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor);
            if(do_link_h_sp_LH_U) data.h_LH_sp.shallowCopy(data.h_U_sp);
            if(do_link_h_sp_UH_L) data.h_UH_sp.shallowCopy(data.h_L_sp);
            data.h_Bperm_sp.resize(data.n_dofs_interface, data.n_dofs_domain, feti.B1[di].nnz);
            if(config->apply_scatter_gather_where == DEVICE::CPU) data.h_applyc_x.resize(data.n_dofs_interface);
            if(config->apply_scatter_gather_where == DEVICE::CPU) data.h_applyc_y.resize(data.n_dofs_interface);
            tm_alloc_host.stop();

            // device memory
            tm_alloc_device.start();
            if(do_alloc_d_sp_L)   data.d_L_sp .resize(data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor);
            if(do_alloc_d_sp_U)   data.d_U_sp .resize(data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor);
            if(do_alloc_d_sp_LH)  data.d_LH_sp.resize(data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor);
            if(do_alloc_d_sp_UH)  data.d_UH_sp.resize(data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor);
            if(do_link_d_sp_LH_U) data.d_LH_sp.shallowCopy(data.d_U_sp);
            if(do_link_d_sp_UH_L) data.d_UH_sp.shallowCopy(data.d_L_sp);
            data.d_Bperm_sp.resize(data.n_dofs_interface, data.n_dofs_domain, feti.B1[di].nnz);
            data.d_apply_x.resize(data.n_dofs_interface);
            data.d_apply_y.resize(data.n_dofs_interface);
            data.d_apply_z.resize(data.n_dofs_domain);
            data.d_apply_w.resize(data.n_dofs_domain);
            data.d_applyg_D2C.resize(data.n_dofs_interface);
            if(do_trsm1_sp) data.buffer_sptrsm1_persistent = gpu::mgm::memalloc_device(data.buffersizes_sptrsm1.persistent);
            if(do_trsm2_sp) data.buffer_sptrsm2_persistent = gpu::mgm::memalloc_device(data.buffersizes_sptrsm2.persistent);
            if(do_trsv1_sp) data.buffer_sptrsv1_persistent = gpu::mgm::memalloc_device(data.buffersizes_sptrsv1.persistent);
            if(do_trsv2_sp) data.buffer_sptrsv2_persistent = gpu::mgm::memalloc_device(data.buffersizes_sptrsv2.persistent);
            if(do_mm)      data.buffer_spmm   = gpu::mgm::memalloc_device(data.buffersize_spmm);
            if(data.buffersize_transL2LH > 0) data.buffer_transL2LH = gpu::mgm::memalloc_device(data.buffersize_transL2LH);
            if(data.buffersize_transU2UH > 0) data.buffer_transU2UH = gpu::mgm::memalloc_device(data.buffersize_transU2UH);
            if(data.buffersize_spmv1 > 0) data.buffer_spmv1 = gpu::mgm::memalloc_device(data.buffersize_spmv1);
            if(data.buffersize_spmv2 > 0) data.buffer_spmv2 = gpu::mgm::memalloc_device(data.buffersize_spmv2);
            if(need_F && data.should_allocate_d_F) d_Fs_allocated[data.allocated_F_index].resize(data.n_dofs_interface + 1, data.n_dofs_interface, data.ld_F);
            tm_alloc_device.stop();
        }
        tm_alloc.stop();

        // prepare matrices on host
        tm_Bperm.start();
        {
            Permutation<I> perm;
            perm.resize(data.n_dofs_domain);
            data.solver_Kreg.getPermutation(perm);
            math::permuteColumns(data.h_Bperm_sp, feti.B1[di], perm);
        }
        tm_Bperm.stop();

        // extract symbolic pattern from the factor
        tm_get_factors.start();
        {
            tm_extract.start();
            if(solver_get_L) data.solver_Kreg.getFactorL(data.h_L_sp, true, false);
            if(solver_get_U) data.solver_Kreg.getFactorU(data.h_U_sp, true, false);
            tm_extract.stop();
            tm_trans_cpu.start();
            if(do_conjtrans_L2LH_h) {
                data.transmap_L2LH.resize(data.h_L_sp.nnz);
                math::csrTranspose(data.h_LH_sp, data.h_L_sp, data.transmap_L2LH, math::CsrTransposeStage::Pattern, true);
            }
            if(do_conjtrans_U2UH_h) {
                data.transmap_U2UH.resize(data.h_U_sp.nnz);
                math::csrTranspose(data.h_UH_sp, data.h_U_sp, data.transmap_U2UH, math::CsrTransposeStage::Pattern, true);
            }
            tm_trans_cpu.stop();
        }
        tm_get_factors.stop();

        // copy some matrices to device
        tm_copyin.start();
        {
            if(do_copyin_L)  gpu::mgm::copy_submit(q, data.d_L_sp,  data.h_L_sp,  true, false);
            if(do_copyin_U)  gpu::mgm::copy_submit(q, data.d_U_sp,  data.h_U_sp,  true, false);
            if(do_copyin_LH) gpu::mgm::copy_submit(q, data.d_LH_sp, data.h_LH_sp, true, false);
            if(do_copyin_UH) gpu::mgm::copy_submit(q, data.d_UH_sp, data.h_UH_sp, true, false);
            gpu::mgm::copy_submit(q, data.d_Bperm_sp, data.h_Bperm_sp, true, true);
            if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        }
        tm_copyin.stop();

        // set the pointers inside the descriptors of some matrices
        tm_setpointers.start();
        {
            if(do_descr_sp_L)  gpu::spblas::descr_matrix_csr_link_data(hs, data.descr_L_sp,  data.d_L_sp);
            if(do_descr_sp_LH) gpu::spblas::descr_matrix_csr_link_data(hs, data.descr_LH_sp, data.d_LH_sp);
            if(do_descr_sp_U)  gpu::spblas::descr_matrix_csr_link_data(hs, data.descr_U_sp,  data.d_U_sp);
            if(do_descr_sp_UH) gpu::spblas::descr_matrix_csr_link_data(hs, data.descr_UH_sp, data.d_UH_sp);
            gpu::spblas::descr_matrix_csr_link_data(hs, data.descr_Bperm_sp, data.d_Bperm_sp);
            if(is_implicit) gpu::spblas::descr_vector_dense_link_data(hs, data.descr_xvec, data.d_apply_x);
            if(is_implicit) gpu::spblas::descr_vector_dense_link_data(hs, data.descr_yvec, data.d_apply_y);
            if(is_implicit) gpu::spblas::descr_vector_dense_link_data(hs, data.descr_zvec, data.d_apply_z);
            if(is_implicit) gpu::spblas::descr_vector_dense_link_data(hs, data.descr_wvec, data.d_apply_w);
        }
        tm_setpointers.stop();

        tm_mainloop_inner.stop();
    }
    tm_mainloop_outer.stop();

    if(need_F)
    {
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_stuff & data = domain_data[di];
            data.d_F.shallowCopy(d_Fs_allocated[data.allocated_F_index]);
            data.d_F.nrows = data.n_dofs_interface;
            data.d_F.ncols = data.n_dofs_interface;
            if((data.hermitian_F_fill == 'L') == (order_F == 'R')) {
                data.d_F.vals += data.d_F.get_ld();
            }
            if(is_path_trsm) {
                gpu::spblas::handle & hs = handles_sparse[di % n_queues];
                gpu::spblas::descr_matrix_dense_link_data(hs, data.descr_F, data.d_F);
            }
        }
    }

    // clean up the mess from buggy openmp in clang
    utils::run_dummy_parallel_region();
    
    // some stuff needed for apply
    tm_applystuff.start();
    {
        if(config->apply_scatter_gather_where == DEVICE::GPU) {
            d_applyg_x_cluster.resize(feti.lambdas.size);
            d_applyg_y_cluster.resize(feti.lambdas.size);
        }
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
        for(size_t di = 0; di < n_domains; di++) {
            h_applyg_xs_pointers.vals[di] = domain_data[di].d_apply_x.vals;
            h_applyg_ys_pointers.vals[di] = domain_data[di].d_apply_y.vals;
            h_applyg_n_dofs_interfaces.vals[di] = domain_data[di].n_dofs_interface;
            h_applyg_D2Cs_pointers.vals[di] = domain_data[di].d_applyg_D2C.vals;
        }
        gpu::mgm::copy_submit(main_q, d_applyg_xs_pointers,       h_applyg_xs_pointers);
        gpu::mgm::copy_submit(main_q, d_applyg_ys_pointers,       h_applyg_ys_pointers);
        gpu::mgm::copy_submit(main_q, d_applyg_n_dofs_interfaces, h_applyg_n_dofs_interfaces);
        gpu::mgm::copy_submit(main_q, d_applyg_D2Cs_pointers,     h_applyg_D2Cs_pointers);
        for(size_t di = 0; di < n_domains; di++) {
            gpu::mgm::copy_submit(main_q, domain_data[di].d_applyg_D2C.vals, feti.D2C[di].data(), feti.D2C[di].size());
        }
        if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
    }
    tm_applystuff.stop();

    // compute mempool requirements
    // to estimate maximum size possibly needed and minimum required
    tm_mempool_reqs.start();
    {
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_stuff & data = domain_data[di];
            size_t & require = data.mempool_requirement_set;
            require = 0;
            require += data.buffersize_tmp_max_set;
        }
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_stuff & data = domain_data[di];
            size_t & require = data.mempool_requirement_update;
            require = 0;
            if(is_explicit && do_alloc_d_dn_L)  require += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
            if(is_explicit && do_alloc_d_dn_U)  require += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
            if(is_explicit && do_alloc_d_dn_LH) require += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
            if(is_explicit && do_alloc_d_dn_UH) require += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
            if(need_X && order_X == 'R') require += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain,    data.n_dofs_interface, data.ld_X);
            if(need_X && order_X == 'C') require += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_interface, data.n_dofs_domain,    data.ld_X);
            if(need_Y && order_X == 'R') require += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain,    data.n_dofs_interface, data.ld_X);
            if(need_Y && order_X == 'C') require += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_interface, data.n_dofs_domain,    data.ld_X);
            if(need_f_tmp) require += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_interface, data.n_dofs_interface, data.ld_F);
            require += data.buffersize_tmp_max_update;
        }
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_stuff & data = domain_data[di];
            size_t & require = data.mempool_requirement_apply;
            require = 0;
            if(is_implicit && do_alloc_d_dn_L)  require += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
            if(is_implicit && do_alloc_d_dn_U)  require += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
            if(is_implicit && do_alloc_d_dn_LH) require += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
            if(is_implicit && do_alloc_d_dn_UH) require += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
            require += data.buffersize_tmp_max_apply;
        }
        mempool_min_size_set = 0;
        mempool_min_size_update = 0;
        mempool_min_size_apply = 0;
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_stuff & data = domain_data[di];
            mempool_min_size_set    = std::max(mempool_min_size_set,    data.mempool_requirement_set);
            mempool_min_size_update = std::max(mempool_min_size_update, data.mempool_requirement_update);
            mempool_min_size_apply  = std::max(mempool_min_size_apply,  data.mempool_requirement_apply);
        }
    }
    tm_mempool_reqs.stop();

    size_t free_mem_before_preprocessing = gpu::mgm::get_device_memory_free();

    // preprocessing
    tm_preprocess.start();
    {
        // smaller mempool just for preprocessing. Some preprocessing functions perform allocations, so I cannot allocate all the memory
        // I allocate at most half of the remaining memory, so the other half is left free. If that is not enough for these allocations, then I don't care, the library should have had better memory management and not allocate by itself
        tm_prepr_poolalloc.start();
        size_t mem_pool_size_request = 0;
        for(size_t qidx = 0; qidx < n_queues; qidx++) {
            size_t mem_pool_size_request_queue = 0;
            for(size_t di = qidx; di < n_domains; di += n_queues) {
                per_domain_stuff & data = domain_data[di];
                size_t mem_pool_size_request_domain = 0;
                mem_pool_size_request_domain += data.buffersize_tmp_max_set;
                mem_pool_size_request_queue = std::max(mem_pool_size_request_queue, mem_pool_size_request_domain);
                // mem_pool_size_request_queue += mem_pool_size_request_domain;
            }
            mem_pool_size_request += mem_pool_size_request_queue;
        }
        mem_pool_size_request = std::min(2 * mem_pool_size_request, gpu::mgm::get_device_memory_free() / 2);
        mem_pool_size_request = std::max(mem_pool_size_request, mempool_min_size_set);
        size_t pool_size_device;
        gpu::mgm::memalloc_device_max(mem_pool_device, pool_size_device, mem_pool_size_request);
        size_t pool_size_device_aligned = (pool_size_device / align_B) * align_B;
        if(pool_size_device_aligned < mempool_min_size_set) {
            eslog::error("not enough memory in mempool for set in rank %d. Need at least %zu MiB, but have only %zu MiB.\n", info::mpi::rank, mempool_min_size_set >> 20, pool_size_device_aligned >> 20);
        }
        cbmba_resource cbmba_res_preprocess(mem_pool_device, pool_size_device_aligned);
        tm_prepr_poolalloc.stop();

        tm_prepr_work.start();
        #pragma omp parallel for schedule(static,1) if(config->concurrency_set == CONCURRENCY::PARALLEL)
        for(size_t di = 0; di < n_domains; di++) {

            gpu::mgm::queue & q = queues[di % n_queues];
            gpu::spblas::handle & hs = handles_sparse[di % n_queues];
            per_domain_stuff & data = domain_data[di];

            gpu::spblas::descr_matrix_dense & descr_X = data.descr_X;
            gpu::spblas::descr_matrix_dense & descr_W = (is_trsm1_outofplace ? data.descr_Y : data.descr_X);
            gpu::spblas::descr_matrix_dense & descr_Z = (is_trsm1_outofplace == is_trsm2_outofplace ? data.descr_X : data.descr_Y);

            void * buffer_tmp = nullptr;

            tm_prepr_allocinpool.start();
            cbmba_res_preprocess.do_transaction([&](){
                buffer_tmp = cbmba_res_preprocess.allocate(data.buffersize_tmp_max_set, align_B);
            });
            tm_prepr_allocinpool.stop();

            // proprocessing stage of the kernels
            tm_prepr_kernels.start();
            {
                tm_trans_gpu.start();
                if(do_conjtrans_L2LH_d) gpu::spblas::transpose<T,I>(hs, data.descr_LH_sp, data.descr_L_sp, true, data.buffersize_transL2LH, data.buffer_transL2LH, 'P');
                if(do_conjtrans_U2UH_d) gpu::spblas::transpose<T,I>(hs, data.descr_UH_sp, data.descr_U_sp, true, data.buffersize_transU2UH, data.buffer_transU2UH, 'P');
                if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                tm_trans_gpu.stop();

                tm_trsm1.start();
                if(do_trsm1_sp_L)  gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_L_sp,  descr_X, descr_W, data.descr_sparse_trsm1, data.buffersizes_sptrsm1, data.buffer_sptrsm1_persistent, buffer_tmp, 'P');
                if(do_trsm1_sp_LH) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_LH_sp, descr_X, descr_W, data.descr_sparse_trsm1, data.buffersizes_sptrsm1, data.buffer_sptrsm1_persistent, buffer_tmp, 'P');
                if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                tm_trsm1.stop();

                if(is_path_trsm) {
                    tm_trsm2.start();
                    if(do_trsm2_sp_U)  gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_U_sp,  descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersizes_sptrsm2, data.buffer_sptrsm2_persistent, buffer_tmp, 'P');
                    if(do_trsm2_sp_UH) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_UH_sp, descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersizes_sptrsm2, data.buffer_sptrsm2_persistent, buffer_tmp, 'P');
                    if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                    tm_trsm2.stop();
                    
                    tm_gemm.start();
                    if(do_mm &&  need_f_tmp) gpu::spblas::mm<T,I>(hs, 'N', 'N', data.descr_Bperm_sp, descr_Z, data.descr_F_tmp, data.buffersize_spmm, data.buffer_spmm, 'P');
                    if(do_mm && !need_f_tmp) gpu::spblas::mm<T,I>(hs, 'N', 'N', data.descr_Bperm_sp, descr_Z, data.descr_F,     data.buffersize_spmm, data.buffer_spmm, 'P');
                    if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                    tm_gemm.stop();
                }

                if(is_implicit)
                {
                    tm_spmv1.start();
                    gpu::spblas::mv<T,I>(hs, 'T', data.descr_Bperm_sp, data.descr_xvec, data.descr_wvec, data.descr_sparse_spmv1, data.buffersize_spmv1, data.buffer_spmv1, 'P');
                    if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                    tm_spmv1.stop();

                    tm_trsv1.start();
                    if(do_trsv1_sp_L)  gpu::spblas::trsv<T,I>(hs, 'N', data.descr_L_sp,  data.descr_wvec, data.descr_zvec, data.descr_sparse_trsv1, data.buffersizes_sptrsv1, data.buffer_sptrsv1_persistent, buffer_tmp, 'P');
                    if(do_trsv1_sp_LH) gpu::spblas::trsv<T,I>(hs, 'H', data.descr_LH_sp, data.descr_wvec, data.descr_zvec, data.descr_sparse_trsv1, data.buffersizes_sptrsv1, data.buffer_sptrsv1_persistent, buffer_tmp, 'P');
                    if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                    tm_trsv1.stop();

                    tm_trsv2.start();
                    if(do_trsv2_sp_U)  gpu::spblas::trsv<T,I>(hs, 'N', data.descr_U_sp,  data.descr_zvec, data.descr_wvec, data.descr_sparse_trsv2, data.buffersizes_sptrsv2, data.buffer_sptrsv2_persistent, buffer_tmp, 'P');
                    if(do_trsv2_sp_UH) gpu::spblas::trsv<T,I>(hs, 'H', data.descr_UH_sp, data.descr_zvec, data.descr_wvec, data.descr_sparse_trsv2, data.buffersizes_sptrsv2, data.buffer_sptrsv2_persistent, buffer_tmp, 'P');
                    if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                    tm_trsv2.stop();

                    tm_spmv2.start();
                    gpu::spblas::mv<T,I>(hs, 'N', data.descr_Bperm_sp, data.descr_wvec, data.descr_yvec, data.descr_sparse_spmv2, data.buffersize_spmv2, data.buffer_spmv2, 'P');
                    if(config->concurrency_set == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                    tm_spmv2.stop();
                }
            }
            tm_prepr_kernels.stop();

            gpu::mgm::submit_host_function(q, [&,buffer_tmp](){
                void * buffer_tmp_ = buffer_tmp;
                cbmba_res_preprocess.deallocate(buffer_tmp_);
            });
        }
        gpu::mgm::device_wait();
        tm_prepr_work.stop();

        tm_prepr_pooldelete.start();
        gpu::mgm::memfree_device(mem_pool_device);
        gpu::mgm::device_wait();
        tm_prepr_pooldelete.stop();
    }
    tm_preprocess.stop();

    size_t free_memory_before_pool = gpu::mgm::get_device_memory_free();
    ssize_t memory_allocated_by_preprocess = (ssize_t)free_mem_before_preprocessing - free_memory_before_pool;

    // memory pool alloc
    tm_poolalloc.start();
    {
        // determine the mempool size such that no blocking would occur
        size_t mem_pool_size_request_update = 0;
        for(size_t qidx = 0; qidx < n_queues; qidx++) {
            size_t mem_pool_size_request_queue = 0;
            for(size_t di = qidx; di < n_domains; di += n_queues) {
                // mem_pool_size_request_queue = std::max(mem_pool_size_request_queue, domain_data[di].mempool_requirement_update);
                mem_pool_size_request_queue += domain_data[di].mempool_requirement_update;
            }
            mem_pool_size_request_update += mem_pool_size_request_queue;
        }
        size_t mem_pool_size_request_apply = 0;
        for(size_t qidx = 0; qidx < n_queues; qidx++) {
            size_t mem_pool_size_request_queue = 0;
            for(size_t di = qidx; di < n_domains; di += n_queues) {
                // mem_pool_size_request_queue = std::max(mem_pool_size_request_queue, domain_data[di].mempool_requirement_apply);
                mem_pool_size_request_queue += domain_data[di].mempool_requirement_apply;
            }
            mem_pool_size_request_apply += mem_pool_size_request_queue;
        }
        size_t mem_pool_size_request = std::max(mem_pool_size_request_update, mem_pool_size_request_apply);
        mem_pool_size_request = 2 * mem_pool_size_request + 1024; // safety margin
        size_t pool_size_device;
        gpu::mgm::memalloc_device_max(mem_pool_device, pool_size_device, mem_pool_size_request);
        // eslog::info("Mem pool requesting %zu MiB, allocated %zu MiB   (%zu B, %zu B)\n", mem_pool_size_request >> 20, pool_size_device >> 20, mem_pool_size_request, pool_size_device);
        size_t pool_size_device_aligned = (pool_size_device / align_B) * align_B;
        cbmba_res_device = std::make_unique<cbmba_resource>(mem_pool_device, pool_size_device_aligned);
    }
    tm_poolalloc.stop();

    if(memory_info_basic) {
        size_t total_mem_sp_factors = 0;
        size_t mem_buffer_sptrsm1 = 0;
        size_t mem_buffer_sptrsm2 = 0;
        size_t mem_buffer_sptrsv1 = 0;
        size_t mem_buffer_sptrsv2 = 0;
        size_t mem_buffer_spmm = 0;
        size_t mem_buffer_transL2LH = 0;
        size_t mem_buffer_transU2UH = 0;
        size_t total_mem_Fs = 0;
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_stuff & data = domain_data[di];
            if(do_alloc_d_sp_L)  total_mem_sp_factors += Matrix_CSR<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor);
            if(do_alloc_d_sp_U)  total_mem_sp_factors += Matrix_CSR<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor);
            if(do_alloc_d_sp_LH) total_mem_sp_factors += Matrix_CSR<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor);
            if(do_alloc_d_sp_UH) total_mem_sp_factors += Matrix_CSR<T,I>::memoryRequirement(data.n_dofs_domain, data.n_dofs_domain, data.n_nz_factor);
            mem_buffer_sptrsm1 += data.buffersizes_sptrsm1.persistent;
            mem_buffer_sptrsm2 += data.buffersizes_sptrsm2.persistent;
            mem_buffer_sptrsv1 += data.buffersizes_sptrsv1.persistent;
            mem_buffer_sptrsv2 += data.buffersizes_sptrsv2.persistent;
            if(is_path_trsm) mem_buffer_spmm += data.buffersize_spmm;
            if(data.buffersize_transL2LH > 0) mem_buffer_transL2LH += data.buffersize_transL2LH;
            if(data.buffersize_transU2UH > 0) mem_buffer_transU2UH += data.buffersize_transU2UH;
            if(data.should_allocate_d_F) total_mem_Fs += Matrix_Dense<T,I>::memoryRequirement(data.n_dofs_interface + 1, data.n_dofs_interface, data.ld_F);
        }
        size_t total_mem_buffers = mem_buffer_sptrsm1 + mem_buffer_sptrsm2 + mem_buffer_sptrsv1 + mem_buffer_sptrsv2 + mem_buffer_spmm + mem_buffer_transL2LH + mem_buffer_transU2UH;

        eslog::info("rank %4d GPU memory capacity [MiB]:                     %9zu\n", info::mpi::rank, gpu::mgm::get_device_memory_capacity() >> 20);
        eslog::info("rank %4d GPU memory used for F matrices [MiB]:          %9zu\n", info::mpi::rank, total_mem_Fs >> 20);
        eslog::info("rank %4d GPU memory used for sparse factors [MiB]:      %9zu\n", info::mpi::rank, total_mem_sp_factors >> 20);
        eslog::info("rank %4d GPU memory used for persistent buffers [MiB]:  %9zu\n", info::mpi::rank, total_mem_buffers >> 20);
        eslog::info("rank %4d   GPU memory used for buffers sptrsm1 [MiB]:   %9zu\n", info::mpi::rank, mem_buffer_sptrsm1 >> 20);
        eslog::info("rank %4d   GPU memory used for buffers sptrsm2 [MiB]:   %9zu\n", info::mpi::rank, mem_buffer_sptrsm2 >> 20);
        eslog::info("rank %4d   GPU memory used for buffers sptrsv1 [MiB]:   %9zu\n", info::mpi::rank, mem_buffer_sptrsv1 >> 20);
        eslog::info("rank %4d   GPU memory used for buffers sptrsv2 [MiB]:   %9zu\n", info::mpi::rank, mem_buffer_sptrsv2 >> 20);
        eslog::info("rank %4d   GPU memory used for buffers spmm [MiB]:      %9zu\n", info::mpi::rank, mem_buffer_spmm >> 20);
        eslog::info("rank %4d   GPU memory used for buffers transL2LH [MiB]: %9zu\n", info::mpi::rank, mem_buffer_transL2LH >> 20);
        eslog::info("rank %4d   GPU memory used for buffers transU2UH [MiB]: %9zu\n", info::mpi::rank, mem_buffer_transU2UH >> 20);
        eslog::info("rank %4d GPU memory allocated in preprocessing [MiB]:   %9zd\n", info::mpi::rank, memory_allocated_by_preprocess >> 20);
        eslog::info("rank %4d GPU memory free for pool and more [MiB]:       %9zu\n", info::mpi::rank, free_memory_before_pool >> 20);
        eslog::info("rank %4d GPU memory used for memory pool [MiB]:         %9zu\n", info::mpi::rank, cbmba_res_device->get_max_capacity() >> 20);
        if(memory_info_detailed) {
            eslog::info("rank %4d   GPU memory need in mempool for set [MiB]:    %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.mempool_requirement_set    >> 20; }).to_string_plain().c_str());
            eslog::info("rank %4d   GPU memory need in mempool for update [MiB]: %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.mempool_requirement_update >> 20; }).to_string_plain().c_str());
            eslog::info("rank %4d   GPU memory need in mempool for apply [MiB]:  %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.mempool_requirement_apply  >> 20; }).to_string_plain().c_str());
            eslog::info("rank %4d GPU memory needed for tmp buffers:\n", info::mpi::rank);
            eslog::info("rank %4d   GPU memory tmp for sptrsm1 preprocess [MiB]: %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.buffersizes_sptrsm1.tmp_preprocess >> 20; }).to_string_plain().c_str());
            eslog::info("rank %4d   GPU memory tmp for sptrsm1 update [MiB]:     %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.buffersizes_sptrsm1.tmp_update     >> 20; }).to_string_plain().c_str());
            eslog::info("rank %4d   GPU memory tmp for sptrsm1 compute [MiB]:    %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.buffersizes_sptrsm1.tmp_compute    >> 20; }).to_string_plain().c_str());
            eslog::info("rank %4d   GPU memory tmp for sptrsm2 preprocess [MiB]: %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.buffersizes_sptrsm2.tmp_preprocess >> 20; }).to_string_plain().c_str());
            eslog::info("rank %4d   GPU memory tmp for sptrsm2 update [MiB]:     %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.buffersizes_sptrsm2.tmp_update     >> 20; }).to_string_plain().c_str());
            eslog::info("rank %4d   GPU memory tmp for sptrsm2 compute [MiB]:    %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.buffersizes_sptrsm2.tmp_compute    >> 20; }).to_string_plain().c_str());
            eslog::info("rank %4d   GPU memory tmp for sptrsm1 preprocess [MiB]: %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.buffersizes_sptrsv1.tmp_preprocess >> 20; }).to_string_plain().c_str());
            eslog::info("rank %4d   GPU memory tmp for sptrsv1 update [MiB]:     %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.buffersizes_sptrsv1.tmp_update     >> 20; }).to_string_plain().c_str());
            eslog::info("rank %4d   GPU memory tmp for sptrsv1 compute [MiB]:    %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.buffersizes_sptrsv1.tmp_compute    >> 20; }).to_string_plain().c_str());
            eslog::info("rank %4d   GPU memory tmp for sptrsv2 preprocess [MiB]: %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.buffersizes_sptrsv2.tmp_preprocess >> 20; }).to_string_plain().c_str());
            eslog::info("rank %4d   GPU memory tmp for sptrsv2 update [MiB]:     %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.buffersizes_sptrsv2.tmp_update     >> 20; }).to_string_plain().c_str());
            eslog::info("rank %4d   GPU memory tmp for sptrsv2 compute [MiB]:    %s\n", info::mpi::rank, minmaxavg<size_t>::compute_from_my_rank_only(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.buffersizes_sptrsv2.tmp_compute    >> 20; }).to_string_plain().c_str());
        }
    }

    if(cbmba_res_device->get_max_capacity() < mempool_min_size_update) {
        eslog::error("not enough memory in mempool for update in rank %d. Need at least %zu MiB, but have only %zu MiB.\n", info::mpi::rank, mempool_min_size_update >> 20, cbmba_res_device->get_max_capacity() >> 20);
    }
    if(cbmba_res_device->get_max_capacity() < mempool_min_size_apply) {
        eslog::error("not enough memory in mempool for apply in rank %d. Need at least %zu MiB, but have only %zu MiB.\n", info::mpi::rank, mempool_min_size_apply >> 20, cbmba_res_device->get_max_capacity() >> 20);
    }

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
    print_timer("Set       libinit", tm_libinit);
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
    print_timer("Set           Bperm", tm_Bperm);
    print_timer("Set           get_factors", tm_get_factors);
    print_timer("Set             extract", tm_extract);
    print_timer("Set             trans_cpu", tm_trans_cpu);
    print_timer("Set           copyin", tm_copyin);
    print_timer("Set           setpointers", tm_setpointers);
    print_timer("Set       applystuff", tm_applystuff);
    print_timer("Set       mempool_reqs", tm_mempool_reqs);
    print_timer("Set       preprocess", tm_preprocess);
    print_timer("Set         prepr_poolalloc", tm_prepr_poolalloc);
    print_timer("Set         prepr_work", tm_prepr_work);
    print_timer("Set           prepr_allocinpool", tm_prepr_allocinpool);
    print_timer("Set           prepr_kernels", tm_prepr_kernels);
    print_timer("Set             trans_gpu", tm_trans_gpu);
    print_timer("Set             trsm1", tm_trsm1);
    print_timer("Set             trsm2", tm_trsm2);
    print_timer("Set             gemm", tm_gemm);
    print_timer("Set             spmv1", tm_spmv1);
    print_timer("Set             trsv1", tm_trsv1);
    print_timer("Set             trsv2", tm_trsv2);
    print_timer("Set             spmv2", tm_spmv2);
    print_timer("Set         prepr_pooldelete", tm_prepr_pooldelete);
    print_timer("Set       poolalloc", tm_poolalloc);
    print_timer("Set       wait", tm_wait);
}

template <typename T, typename I>
void TotalFETIGpu<T,I>::update(const step::Step &step)
{
    if(stage != 2 && stage != 3) eslog::error("update: invalud order of operations in dualop\n");

    my_timer tm_total(timers_basic), tm_mainloop_outer(timers_basic), tm_compute_d(timers_basic), tm_wait(timers_basic);
    my_timer tm_mainloop_inner(timers_detailed), tm_Kreg_combine(timers_detailed), tm_solver_commit(timers_detailed), tm_fact_numeric(timers_detailed), tm_get_factors(timers_detailed), tm_extract(timers_detailed), tm_trans_cpu(timers_detailed), tm_allocinpool(timers_detailed), tm_setpointers(timers_detailed), tm_copyin(timers_detailed), tm_trans_gpu(timers_detailed), tm_descr_update(timers_detailed), tm_descr_update_trsm1(timers_detailed), tm_descr_update_trsm2(timers_detailed), tm_descr_update_trsv1(timers_detailed), tm_descr_update_trsv2(timers_detailed), tm_sp2dn(timers_detailed), tm_kernels_compute(timers_detailed), tm_trsm1(timers_detailed), tm_trsm2(timers_detailed), tm_gemm(timers_detailed), tm_fcopy(timers_detailed), tm_syrk(timers_detailed), tm_freeinpool(timers_detailed), tm_freeinpool_exec(timers_detailed);

    tm_total.start();

    gpu::mgm::set_device(device);

    tm_mainloop_outer.start();
    #pragma omp parallel for schedule(static,1) if(config->concurrency_update == CONCURRENCY::PARALLEL)
    for(size_t di = 0; di < n_domains; di++) {
        tm_mainloop_inner.start();

        gpu::mgm::queue & q = queues[di % n_queues];
        gpu::dnblas::handle & hd = handles_dense[di % n_queues];
        gpu::spblas::handle & hs = handles_sparse[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        std::unique_ptr<Matrix_Dense<T,I,cbmba_d>> & d_X = data.d_X;
        std::unique_ptr<Matrix_Dense<T,I,cbmba_d>> & d_W = (is_trsm1_outofplace ? data.d_Y : data.d_X);
        // std::unique_ptr<Matrix_Dense<T,I,cbmba_d>> & d_Z = (is_trsm1_outofplace == is_trsm2_outofplace ? data.d_X : data.d_Y);
        gpu::spblas::descr_matrix_dense & descr_X = data.descr_X;
        gpu::spblas::descr_matrix_dense & descr_W = (is_trsm1_outofplace ? data.descr_Y : data.descr_X);
        gpu::spblas::descr_matrix_dense & descr_Z = (is_trsm1_outofplace == is_trsm2_outofplace ? data.descr_X : data.descr_Y);

        void * buffer_tmp = nullptr;

        // Kreg = K + RegMat numeric values
        tm_Kreg_combine.start();
        {
            math::sumCombined(data.Kreg, T{1.0}, feti.K[di], feti.RegMat[di]);
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
            tm_trans_cpu.start();
            if(do_conjtrans_L2LH_h) math::csrTranspose(data.h_LH_sp, data.h_L_sp, data.transmap_L2LH, math::CsrTransposeStage::Values, true);
            if(do_conjtrans_U2UH_h) math::csrTranspose(data.h_UH_sp, data.h_U_sp, data.transmap_U2UH, math::CsrTransposeStage::Values, true);
            tm_trans_cpu.stop();
        }
        tm_get_factors.stop();

        // temporary allocations using the memory pool
        tm_allocinpool.start();
        {
            if(is_explicit) {
                cbmba_d ator_d(*cbmba_res_device, align_B);
                cbmba_res_device->do_transaction([&](){
                    data.d_L_dn  = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
                    data.d_U_dn  = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
                    data.d_LH_dn = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
                    data.d_UH_dn = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
                    data.d_X = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
                    data.d_Y = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
                    data.d_F_tmp = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);

                    if(do_alloc_d_dn_L)   data.d_L_dn ->resize(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
                    if(do_alloc_d_dn_U)   data.d_U_dn ->resize(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
                    if(do_alloc_d_dn_LH)  data.d_LH_dn->resize(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
                    if(do_alloc_d_dn_UH)  data.d_UH_dn->resize(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
                    if(do_link_d_dn_LH_U) data.d_LH_dn->shallowCopy(*data.d_U_dn);
                    if(do_link_d_dn_UH_L) data.d_UH_dn->shallowCopy(*data.d_L_dn);
                    if(order_X == 'R')           data.d_X->resize(data.n_dofs_domain,    data.n_dofs_interface, data.ld_X);
                    if(order_X == 'C')           data.d_X->resize(data.n_dofs_interface, data.n_dofs_domain,    data.ld_X);
                    if(order_X == 'R' && need_Y) data.d_Y->resize(data.n_dofs_domain,    data.n_dofs_interface, data.ld_X);
                    if(order_X == 'C' && need_Y) data.d_Y->resize(data.n_dofs_interface, data.n_dofs_domain,    data.ld_X);
                    if(need_f_tmp) data.d_F_tmp->resize(data.n_dofs_interface, data.n_dofs_interface, data.ld_F);
                    buffer_tmp = cbmba_res_device->allocate(data.buffersize_tmp_max_update, align_B);
                });
            }
        }
        tm_allocinpool.stop();

        if(is_explicit) {
            gpu::dnblas::buffer_set(hd, buffer_tmp, data.buffersize_tmp_max_update);
        }

        // copy the new factors to device
        tm_copyin.start();
        {
            if(do_copyin_L)  gpu::mgm::copy_submit(q, data.d_L_sp,  data.h_L_sp,  false, true);
            if(do_copyin_U)  gpu::mgm::copy_submit(q, data.d_U_sp,  data.h_U_sp,  false, true);
            if(do_copyin_LH) gpu::mgm::copy_submit(q, data.d_LH_sp, data.h_LH_sp, false, true);
            if(do_copyin_UH) gpu::mgm::copy_submit(q, data.d_UH_sp, data.h_UH_sp, false, true);
            if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        }
        tm_copyin.stop();

        // set the pointers inside the descriptors of the rest of the matrices
        tm_setpointers.start();
        {
            if(is_explicit) {
                if(do_descr_dn_L)  gpu::spblas::descr_matrix_dense_link_data(hs, data.descr_L_dn,  *data.d_L_dn);
                if(do_descr_dn_LH) gpu::spblas::descr_matrix_dense_link_data(hs, data.descr_LH_dn, *data.d_LH_dn);
                if(do_descr_dn_U)  gpu::spblas::descr_matrix_dense_link_data(hs, data.descr_U_dn,  *data.d_U_dn);
                if(do_descr_dn_UH) gpu::spblas::descr_matrix_dense_link_data(hs, data.descr_UH_dn, *data.d_UH_dn);
                if(need_X) gpu::spblas::descr_matrix_dense_link_data(hs, data.descr_X, *data.d_X);
                if(need_Y) gpu::spblas::descr_matrix_dense_link_data(hs, data.descr_Y, *data.d_Y);
                if(need_f_tmp) gpu::spblas::descr_matrix_dense_link_data(hs, data.descr_F_tmp, *data.d_F_tmp);
            }
            if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        }
        tm_setpointers.stop();

        // transpose csr matrices on device
        tm_trans_gpu.start();
        {
            if(do_conjtrans_L2LH_d) gpu::spblas::transpose<T,I>(hs, data.descr_LH_sp, data.descr_L_sp, true, data.buffersize_transL2LH, data.buffer_transL2LH, 'C');
            if(do_conjtrans_U2UH_d) gpu::spblas::transpose<T,I>(hs, data.descr_UH_sp, data.descr_U_sp, true, data.buffersize_transU2UH, data.buffer_transU2UH, 'C');
            if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        }
        tm_trans_gpu.stop();

        // update sparse trsm descriptors to reflect the new matrix values, possibly re-preprocess
        tm_descr_update.start();
        {
            tm_descr_update_trsm1.start();
            if(do_trsm1_sp_L)  gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_L_sp,  descr_X, descr_W, data.descr_sparse_trsm1, data.buffersizes_sptrsm1, data.buffer_sptrsm1_persistent, buffer_tmp, 'U');
            if(do_trsm1_sp_LH) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_LH_sp, descr_X, descr_W, data.descr_sparse_trsm1, data.buffersizes_sptrsm1, data.buffer_sptrsm1_persistent, buffer_tmp, 'U');
            if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
            tm_descr_update_trsm1.stop();
            tm_descr_update_trsm2.start();
            if(do_trsm2_sp_U)  gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_U_sp,  descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersizes_sptrsm2, data.buffer_sptrsm2_persistent, buffer_tmp, 'U');
            if(do_trsm2_sp_UH) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_UH_sp, descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersizes_sptrsm2, data.buffer_sptrsm2_persistent, buffer_tmp, 'U');
            if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
            tm_descr_update_trsm2.stop();
            tm_descr_update_trsv1.start();
            if(do_trsv1_sp_L)  gpu::spblas::trsv<T,I>(hs, 'N', data.descr_L_sp,  data.descr_wvec, data.descr_zvec, data.descr_sparse_trsv1, data.buffersizes_sptrsv1, data.buffer_sptrsv1_persistent, buffer_tmp, 'U');
            if(do_trsv1_sp_LH) gpu::spblas::trsv<T,I>(hs, 'H', data.descr_LH_sp, data.descr_wvec, data.descr_zvec, data.descr_sparse_trsv1, data.buffersizes_sptrsv1, data.buffer_sptrsv1_persistent, buffer_tmp, 'U');
            if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
            tm_descr_update_trsv1.stop();
            tm_descr_update_trsv2.start();
            if(do_trsv2_sp_U)  gpu::spblas::trsv<T,I>(hs, 'N', data.descr_U_sp,  data.descr_zvec, data.descr_wvec, data.descr_sparse_trsv2, data.buffersizes_sptrsv2, data.buffer_sptrsv2_persistent, buffer_tmp, 'U');
            if(do_trsv2_sp_UH) gpu::spblas::trsv<T,I>(hs, 'H', data.descr_UH_sp, data.descr_zvec, data.descr_wvec, data.descr_sparse_trsv2, data.buffersizes_sptrsv2, data.buffer_sptrsv2_persistent, buffer_tmp, 'U');
            if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
            tm_descr_update_trsv2.stop();
        }
        tm_descr_update.stop();

        // sparse to dense on device
        tm_sp2dn.start();
        {
            if(is_explicit) {
                if(do_sp2dn_L)  gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_L_sp,  data.descr_L_dn,  data.buffersize_tmp_max_update, buffer_tmp, 'C');
                if(do_sp2dn_U)  gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_U_sp,  data.descr_U_dn,  data.buffersize_tmp_max_update, buffer_tmp, 'C');
                if(do_sp2dn_LH) gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_LH_sp, data.descr_LH_dn, data.buffersize_tmp_max_update, buffer_tmp, 'C');
                if(do_sp2dn_UH) gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_UH_sp, data.descr_UH_dn, data.buffersize_tmp_max_update, buffer_tmp, 'C');
                if(need_X) gpu::spblas::sparse_to_dense<T,I>(hs, 'T', data.descr_Bperm_sp, data.descr_X, data.buffersize_tmp_max_update, buffer_tmp, 'C');
                if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
            }
        }
        tm_sp2dn.stop();

        if(is_explicit) {
            // perform the actual assembly
            tm_kernels_compute.start();
            {
                tm_trsm1.start();
                if(do_trsm1_sp_L)  gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_L_sp,  descr_X, descr_W, data.descr_sparse_trsm1, data.buffersizes_sptrsm1, data.buffer_sptrsm1_persistent, buffer_tmp, 'C');
                if(do_trsm1_sp_LH) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_LH_sp, descr_X, descr_W, data.descr_sparse_trsm1, data.buffersizes_sptrsm1, data.buffer_sptrsm1_persistent, buffer_tmp, 'C');
                if(do_trsm1_dn_L)  gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, data.d_L_dn->vals,  data.ld_domain, 'R', 'N', 'L', d_X->vals, data.ld_X, order_X, 'N');
                if(do_trsm1_dn_LH) gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, data.d_LH_dn->vals, data.ld_domain, 'R', 'H', 'U', d_X->vals, data.ld_X, order_X, 'N');
                if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                tm_trsm1.stop();

                if(is_path_herk) {
                    tm_syrk.start();
                    gpu::dnblas::herk<T,I>(hd, data.n_dofs_interface, data.n_dofs_domain, d_W->vals, data.ld_X, order_X, 'H', data.d_F.vals, data.ld_F, order_F, data.hermitian_F_fill);
                    if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                    tm_syrk.stop();
                }
                if(is_path_trsm) {
                    tm_trsm2.start();
                    if(do_trsm2_sp_U)  gpu::spblas::trsm<T,I>(hs, 'N', 'N', 'N', data.descr_U_sp,  descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersizes_sptrsm2, data.buffer_sptrsm2_persistent, buffer_tmp, 'C');
                    if(do_trsm2_sp_UH) gpu::spblas::trsm<T,I>(hs, 'H', 'N', 'N', data.descr_UH_sp, descr_W, descr_Z, data.descr_sparse_trsm2, data.buffersizes_sptrsm2, data.buffer_sptrsm2_persistent, buffer_tmp, 'C');
                    if(do_trsm2_dn_U)  gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, data.d_U_dn->vals,  data.ld_domain, 'R', 'N', 'U', d_W->vals, data.ld_X, order_X, 'N');
                    if(do_trsm2_dn_UH) gpu::dnblas::trsm<T,I>(hd, 'L', data.n_dofs_domain, data.n_dofs_interface, data.d_UH_dn->vals, data.ld_domain, 'R', 'H', 'L', d_W->vals, data.ld_X, order_X, 'N');
                    if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                    tm_trsm2.stop();

                    tm_gemm.start();
                    if( need_f_tmp) gpu::spblas::mm<T,I>(hs, 'N', 'N', data.descr_Bperm_sp, descr_Z, data.descr_F_tmp, data.buffersize_spmm, data.buffer_spmm, 'C');
                    if(!need_f_tmp) gpu::spblas::mm<T,I>(hs, 'N', 'N', data.descr_Bperm_sp, descr_Z, data.descr_F,     data.buffersize_spmm, data.buffer_spmm, 'C');
                    if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                    tm_gemm.stop();

                    tm_fcopy.start();
                    if(need_f_tmp) gpu::kernels::copy_matrix_triangle(q, data.d_F, *data.d_F_tmp, data.hermitian_F_fill, order_F);
                    if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
                    tm_fcopy.stop();
                }
            }
            tm_kernels_compute.stop();
        }

        gpu::dnblas::buffer_unset(hd);

        // free the temporary memory from the pool
        tm_freeinpool.start();
        if(is_explicit)
        {
            gpu::mgm::submit_host_function(q, [&,di,buffer_tmp]() {
                domain_data[di].d_L_dn.reset();
                domain_data[di].d_LH_dn.reset();
                domain_data[di].d_U_dn.reset();
                domain_data[di].d_UH_dn.reset();
                domain_data[di].d_X.reset();
                domain_data[di].d_Y.reset();
                domain_data[di].d_F_tmp.reset();
                void * buffer_tmp_ = buffer_tmp;
                cbmba_res_device->deallocate(buffer_tmp_);
            });
            if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        }
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
        for(size_t di = 0; di < n_domains; di++) {
            domain_data[di].solver_Kreg.solve(feti.f[di], Kplus_fs[di]);
        }
        applyB(feti, Kplus_fs, d);
        d.synchronize();
        math::add(d, T{-1}, feti.c);
    }
    tm_compute_d.stop();

    tm_wait.start();
    if(config->synchronize_update) gpu::mgm::device_wait();
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
    print_timer("Update          trans_cpu", tm_trans_cpu);
    print_timer("Update        allocinpool", tm_allocinpool);
    print_timer("Update        copyin", tm_copyin);
    print_timer("Update        setpointers", tm_setpointers);
    print_timer("Update        trans_gpu", tm_trans_gpu);
    print_timer("Update        descr_update", tm_descr_update);
    print_timer("Update          descr_update_trsm1", tm_descr_update_trsm1);
    print_timer("Update          descr_update_trsm2", tm_descr_update_trsm2);
    print_timer("Update          descr_update_trsv1", tm_descr_update_trsv1);
    print_timer("Update          descr_update_trsv2", tm_descr_update_trsv2);
    print_timer("Update        sp2dn", tm_sp2dn);
    print_timer("Update        kernels_compute", tm_kernels_compute);
    print_timer("Update          trsm1", tm_trsm1);
    print_timer("Update          trsm2", tm_trsm2);
    print_timer("Update          gemm", tm_gemm);
    print_timer("Update          fcopy", tm_fcopy);
    print_timer("Update          syrk", tm_syrk);
    print_timer("Update        freeinpool", tm_freeinpool);
    print_timer("Update    compute_d", tm_compute_d);
    print_timer("Update    wait", tm_wait);

    print(step);
}

template<typename T, typename I>
void TotalFETIGpu<T,I>::apply_explicit_sgcpu(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    my_timer tm_total(timers_basic);
    my_timer tm_bufferalloc(timers_detailed), tm_apply_outer(timers_detailed), tm_apply_inner(timers_detailed), tm_scatter(timers_detailed), tm_copyin(timers_detailed), tm_mv(timers_detailed), tm_copyout(timers_detailed), tm_bufferfree(timers_detailed), tm_zerofill(timers_detailed), tm_wait(timers_detailed), tm_gather(timers_detailed), tm_gather_inner(timers_detailed);

    tm_total.start();

    tm_bufferalloc.start();
    std::vector<void*> buffers_tmp(n_queues);
    for(size_t qi = 0; qi < n_queues; qi++) {
        size_t tmp_buffer_needed_queue = 0;
        for(size_t di = qi; di < n_domains; di += n_queues) {
            tmp_buffer_needed_queue = std::max(tmp_buffer_needed_queue, domain_data[di].buffersize_tmp_max_apply);
        }
        buffers_tmp[qi] = cbmba_res_device->allocate(tmp_buffer_needed_queue, align_B);
        gpu::dnblas::buffer_set(handles_dense[qi], buffers_tmp[qi], tmp_buffer_needed_queue);
    }
    tm_bufferalloc.stop();

    tm_apply_outer.start();
    #pragma omp parallel for schedule(static,1) if(config->concurrency_apply == CONCURRENCY::PARALLEL)
    for(size_t di = 0; di < n_domains; di++) {
        tm_apply_inner.start();

        gpu::mgm::queue & q = queues[di % n_queues];
        gpu::dnblas::handle & hd = handles_dense[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        tm_scatter.start();
        for(I i = 0; i < data.n_dofs_interface; i++) {
            data.h_applyc_x.vals[i] = x_cluster.vals[feti.D2C[di][i]];
        }
        tm_scatter.stop();

        tm_copyin.start();
        gpu::mgm::copy_submit(q, data.d_apply_x, data.h_applyc_x);
        if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        tm_copyin.stop();

        tm_mv.start();
        if( is_system_hermitian) gpu::dnblas::hemv<T,I>(hd, data.d_F.nrows, data.d_F.vals, data.d_F.get_ld(), order_F, 'N', data.hermitian_F_fill, data.d_apply_x.vals, data.d_apply_y.vals);
        if(!is_system_hermitian) gpu::dnblas::gemv<T,I>(hd, data.d_F.nrows, data.d_F.ncols, data.d_F.vals, data.d_F.get_ld(), order_F, 'N', data.d_apply_x.vals, data.d_apply_y.vals);
        if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        tm_mv.stop();

        tm_copyout.start();
        gpu::mgm::copy_submit(q, data.h_applyc_y, data.d_apply_y);
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

    tm_bufferfree.start();
    for(size_t qi = 0; qi < n_queues; qi++) {
        gpu::dnblas::buffer_unset(handles_dense[qi]);
        cbmba_res_device->deallocate(buffers_tmp[qi]);
    }
    tm_bufferfree.stop();

    tm_gather.start();
    #pragma omp parallel for schedule(static,1) if(config->concurrency_apply == CONCURRENCY::PARALLEL)
    for(size_t di = 0; di < n_domains; di++) {
        tm_gather_inner.start();
        per_domain_stuff & data = domain_data[di];
        for(I i = 0; i < data.n_dofs_interface; i++) {
            if constexpr(utils::is_real<T>()) {
                #pragma omp atomic
                y_cluster.vals[feti.D2C[di][i]] += data.h_applyc_y.vals[i];
            }
            if constexpr(utils::is_complex<T>()) {
                #pragma omp atomic
                utils::real_ref(y_cluster.vals[feti.D2C[di][i]]) += utils::real_ref(data.h_applyc_y.vals[i]);
                #pragma omp atomic
                utils::imag_ref(y_cluster.vals[feti.D2C[di][i]]) += utils::imag_ref(data.h_applyc_y.vals[i]);
            }
        }
        tm_gather_inner.stop();
    }
    tm_gather.stop();

    tm_total.stop();

    print_timer("Apply   total", tm_total);
    print_timer("Apply     bufferalloc", tm_bufferalloc);
    print_timer("Apply     apply_outer", tm_apply_outer);
    print_timer("Apply       apply_inner", tm_apply_inner);
    print_timer("Apply         scatter", tm_scatter);
    print_timer("Apply         copyin", tm_copyin);
    print_timer("Apply         mv", tm_mv);
    print_timer("Apply         copyout", tm_copyout);
    print_timer("Apply     zerofill", tm_zerofill);
    print_timer("Apply     wait", tm_wait);
    print_timer("Apply     bufferfree", tm_bufferfree);
    print_timer("Apply     gather", tm_gather);
    print_timer("Apply     gather_inner", tm_gather_inner);
}

template<typename T, typename I>
void TotalFETIGpu<T,I>::apply_explicit_sggpu(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    my_timer tm_total(timers_basic);
    my_timer tm_bufferalloc(timers_detailed), tm_copyin(timers_detailed), tm_scatter(timers_detailed), tm_mv_outer(timers_detailed), tm_mv(timers_detailed), tm_zerofill(timers_detailed), tm_gather(timers_detailed), tm_copyout(timers_detailed), tm_wait(timers_detailed), tm_bufferfree(timers_detailed);

    tm_total.start();

    tm_bufferalloc.start();
    std::vector<void*> buffers_tmp(n_queues);
    for(size_t qi = 0; qi < n_queues; qi++) {
        size_t tmp_buffer_needed_queue = 0;
        for(size_t di = qi; di < n_domains; di += n_queues) {
            tmp_buffer_needed_queue = std::max(tmp_buffer_needed_queue, domain_data[di].buffersize_tmp_max_apply);
        }
        buffers_tmp[qi] = cbmba_res_device->allocate(tmp_buffer_needed_queue, align_B);
        gpu::dnblas::buffer_set(handles_dense[qi], buffers_tmp[qi], tmp_buffer_needed_queue);
    }
    tm_bufferalloc.stop();

    // copy x_cluster to device
    tm_copyin.start();
    gpu::mgm::copy_submit(main_q, d_applyg_x_cluster, x_cluster);
    if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
    tm_copyin.stop();

    // scatter
    tm_scatter.start();
    gpu::kernels::DCmap_scatter(main_q, d_applyg_xs_pointers, d_applyg_n_dofs_interfaces, d_applyg_x_cluster, d_applyg_D2Cs_pointers);
    if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
    tm_scatter.stop();

    gpu::mgm::queue_async_barrier({main_q}, queues);

    // apply
    tm_mv_outer.start();
    #pragma omp parallel for schedule(static,1) if(config->concurrency_apply == CONCURRENCY::PARALLEL)
    for(size_t di = 0; di < n_domains; di++) {
        gpu::mgm::queue & q = queues[di % n_queues];
        gpu::dnblas::handle & hd = handles_dense[di % n_queues];
        per_domain_stuff & data = domain_data[di];

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
    gpu::mgm::copy_submit(main_q, y_cluster, d_applyg_y_cluster);
    if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
    tm_copyout.stop();

    // wait
    tm_wait.start();
    gpu::mgm::device_wait();
    tm_wait.stop();

    tm_bufferfree.start();
    for(size_t qi = 0; qi < n_queues; qi++) {
        gpu::dnblas::buffer_unset(handles_dense[qi]);
        cbmba_res_device->deallocate(buffers_tmp[qi]);
    }
    tm_bufferfree.stop();

    tm_total.stop();

    print_timer("Apply   total", tm_total);
    print_timer("Apply     bufferalloc", tm_bufferalloc);
    print_timer("Apply     copyin", tm_copyin);
    print_timer("Apply     scatter", tm_scatter);
    print_timer("Apply     mv_outer", tm_mv_outer);
    print_timer("Apply       mv", tm_mv);
    print_timer("Apply     zerofill", tm_zerofill);
    print_timer("Apply     gather", tm_gather);
    print_timer("Apply     copyout", tm_copyout);
    print_timer("Apply     wait", tm_wait);
    print_timer("Apply     bufferfree", tm_bufferfree);
}

template<typename T, typename I>
void TotalFETIGpu<T,I>::apply_implicit_sgcpu(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    my_timer tm_total(timers_basic);
    my_timer tm_apply_outer(timers_detailed), tm_apply_inner(timers_detailed), tm_scatter(timers_detailed), tm_copyin(timers_detailed), tm_compute(timers_detailed), tm_copyout(timers_detailed), tm_zerofill(timers_detailed), tm_wait(timers_detailed), tm_gather(timers_detailed), tm_gather_inner(timers_detailed), tm_allocinpool(timers_detailed), tm_freeinpool(timers_detailed), tm_setpointers(timers_detailed), tm_sp2dn(timers_detailed), tm_spmv1(timers_detailed), tm_spmv2(timers_detailed), tm_trsv1(timers_detailed), tm_trsv2(timers_detailed);

    tm_total.start();

    tm_apply_outer.start();
    #pragma omp parallel for schedule(static,1) if(config->concurrency_apply == CONCURRENCY::PARALLEL)
    for(size_t di = 0; di < n_domains; di++) {
        tm_apply_inner.start();

        gpu::mgm::queue & q = queues[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        tm_scatter.start();
        for(I i = 0; i < data.n_dofs_interface; i++) {
            data.h_applyc_x.vals[i] = x_cluster.vals[feti.D2C[di][i]];
        }
        tm_scatter.stop();

        tm_copyin.start();
        gpu::mgm::copy_submit(q, data.d_apply_x, data.h_applyc_x);
        if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        tm_copyin.stop();

        apply_implicit_compute(di, tm_allocinpool, tm_freeinpool, tm_setpointers, tm_sp2dn, tm_compute, tm_spmv1, tm_trsv1, tm_trsv2, tm_spmv2);

        tm_copyout.start();
        gpu::mgm::copy_submit(q, data.h_applyc_y, data.d_apply_y);
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
    for(size_t di = 0; di < n_domains; di++) {
        tm_gather_inner.start();
        per_domain_stuff & data = domain_data[di];
        for(I i = 0; i < data.n_dofs_interface; i++) {
            if constexpr(utils::is_real<T>()) {
                #pragma omp atomic
                y_cluster.vals[feti.D2C[di][i]] += data.h_applyc_y.vals[i];
            }
            if constexpr(utils::is_complex<T>()) {
                #pragma omp atomic
                utils::real_ref(y_cluster.vals[feti.D2C[di][i]]) += utils::real_ref(data.h_applyc_y.vals[i]);
                #pragma omp atomic
                utils::imag_ref(y_cluster.vals[feti.D2C[di][i]]) += utils::imag_ref(data.h_applyc_y.vals[i]);
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
    print_timer("Apply         allocinpool", tm_allocinpool);
    print_timer("Apply         setpointers", tm_setpointers);
    print_timer("Apply         sp2dn", tm_sp2dn);
    print_timer("Apply         compute", tm_compute);
    print_timer("Apply           spmv1", tm_spmv1);
    print_timer("Apply           trsv1", tm_trsv1);
    print_timer("Apply           trsv2", tm_trsv2);
    print_timer("Apply           spmv2", tm_spmv2);
    print_timer("Apply         freeinpool", tm_freeinpool);
    print_timer("Apply         copyout", tm_copyout);
    print_timer("Apply     zerofill", tm_zerofill);
    print_timer("Apply     wait", tm_wait);
    print_timer("Apply     gather", tm_gather);
    print_timer("Apply     gather_inner", tm_gather_inner);
}

template<typename T, typename I>
void TotalFETIGpu<T,I>::apply_implicit_sggpu(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    my_timer tm_total(timers_basic);
    my_timer tm_copyin(timers_detailed), tm_scatter(timers_detailed), tm_apply_outer(timers_detailed), tm_compute(timers_detailed), tm_zerofill(timers_detailed), tm_gather(timers_detailed), tm_copyout(timers_detailed), tm_wait(timers_detailed), tm_allocinpool(timers_detailed), tm_freeinpool(timers_detailed), tm_setpointers(timers_detailed), tm_sp2dn(timers_detailed), tm_spmv1(timers_detailed), tm_spmv2(timers_detailed), tm_trsv1(timers_detailed), tm_trsv2(timers_detailed);

    tm_total.start();

    // copy x_cluster to device
    tm_copyin.start();
    gpu::mgm::copy_submit(main_q, d_applyg_x_cluster, x_cluster);
    if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
    tm_copyin.stop();

    // scatter
    tm_scatter.start();
    gpu::kernels::DCmap_scatter(main_q, d_applyg_xs_pointers, d_applyg_n_dofs_interfaces, d_applyg_x_cluster, d_applyg_D2Cs_pointers);
    if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(main_q);
    tm_scatter.stop();

    // apply
    tm_apply_outer.start();
    gpu::mgm::queue_async_barrier({main_q}, queues);
    #pragma omp parallel for schedule(static,1) if(config->concurrency_apply == CONCURRENCY::PARALLEL)
    for(size_t di = 0; di < n_domains; di++) {
        apply_implicit_compute(di, tm_allocinpool, tm_freeinpool, tm_setpointers, tm_sp2dn, tm_compute, tm_spmv1, tm_trsv1, tm_trsv2, tm_spmv2);
    }
    tm_apply_outer.stop();

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
    gpu::mgm::copy_submit(main_q, y_cluster, d_applyg_y_cluster);
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
    print_timer("Apply     apply_outer", tm_apply_outer);
    print_timer("Apply       allocinpool", tm_allocinpool);
    print_timer("Apply       setpointers", tm_setpointers);
    print_timer("Apply       sp2dn", tm_sp2dn);
    print_timer("Apply       compute", tm_compute);
    print_timer("Apply         spmv1", tm_spmv1);
    print_timer("Apply         trsv1", tm_trsv1);
    print_timer("Apply         trsv2", tm_trsv2);
    print_timer("Apply         spmv2", tm_spmv2);
    print_timer("Apply       freeinpool", tm_freeinpool);
    print_timer("Apply     zerofill", tm_zerofill);
    print_timer("Apply     gather", tm_gather);
    print_timer("Apply     copyout", tm_copyout);
    print_timer("Apply     wait", tm_wait);
}

template<typename T, typename I>
void TotalFETIGpu<T,I>::apply_implicit_compute(size_t di, my_timer & tm_allocinpool, my_timer & tm_freeinpool, my_timer & tm_setpointers, my_timer & tm_sp2dn, my_timer & tm_compute, my_timer & tm_spmv1, my_timer & tm_trsv1, my_timer & tm_trsv2, my_timer & tm_spmv2)
{
    gpu::mgm::queue & q = queues[di % n_queues];
    gpu::dnblas::handle & hd = handles_dense[di % n_queues];
    gpu::spblas::handle & hs = handles_sparse[di % n_queues];
    per_domain_stuff & data = domain_data[di];

    void * buffer_tmp = nullptr;

    tm_allocinpool.start();
    {
        cbmba_d ator_d(*cbmba_res_device, align_B);
        cbmba_res_device->do_transaction([&](){
            data.d_L_dn  = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
            data.d_U_dn  = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
            data.d_LH_dn = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
            data.d_UH_dn = std::make_unique<Matrix_Dense<T,I,cbmba_d>>(ator_d);
            if(do_alloc_d_dn_L)   data.d_L_dn ->resize(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
            if(do_alloc_d_dn_U)   data.d_U_dn ->resize(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
            if(do_alloc_d_dn_LH)  data.d_LH_dn->resize(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
            if(do_alloc_d_dn_UH)  data.d_UH_dn->resize(data.n_dofs_domain, data.n_dofs_domain, data.ld_domain);
            if(do_link_d_dn_LH_U) data.d_LH_dn->shallowCopy(*data.d_U_dn);
            if(do_link_d_dn_UH_L) data.d_UH_dn->shallowCopy(*data.d_L_dn);

            buffer_tmp = cbmba_res_device->allocate(data.buffersize_tmp_max_apply, align_B);
        });
    }
    tm_allocinpool.stop();

    gpu::dnblas::buffer_set(hd, buffer_tmp, data.buffersize_tmp_max_apply);

    tm_setpointers.start();
    {
        if(do_descr_dn_L)  gpu::spblas::descr_matrix_dense_link_data(hs, data.descr_L_dn,  *data.d_L_dn);
        if(do_descr_dn_LH) gpu::spblas::descr_matrix_dense_link_data(hs, data.descr_LH_dn, *data.d_LH_dn);
        if(do_descr_dn_U)  gpu::spblas::descr_matrix_dense_link_data(hs, data.descr_U_dn,  *data.d_U_dn);
        if(do_descr_dn_UH) gpu::spblas::descr_matrix_dense_link_data(hs, data.descr_UH_dn, *data.d_UH_dn);
    }
    tm_setpointers.stop();

    tm_sp2dn.start();
    {
        if(do_sp2dn_L)  gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_L_sp,  data.descr_L_dn,  data.buffersize_tmp_max_apply, buffer_tmp, 'C');
        if(do_sp2dn_U)  gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_U_sp,  data.descr_U_dn,  data.buffersize_tmp_max_apply, buffer_tmp, 'C');
        if(do_sp2dn_LH) gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_LH_sp, data.descr_LH_dn, data.buffersize_tmp_max_apply, buffer_tmp, 'C');
        if(do_sp2dn_UH) gpu::spblas::sparse_to_dense<T,I>(hs, 'N', data.descr_UH_sp, data.descr_UH_dn, data.buffersize_tmp_max_apply, buffer_tmp, 'C');
        if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
    }
    tm_sp2dn.stop();

    tm_compute.start();
    {
        Vector_Dense<T,I,Ad> & vec_1 = data.d_apply_w;
        Vector_Dense<T,I,Ad> & vec_2 = (is_factor1_sparse ? data.d_apply_z : data.d_apply_w);
        // Vector_Dense<T,I,Ad> & vec_3 = ((is_factor1_sparse == is_factor2_sparse) ? data.d_apply_w : data.d_apply_z);
        gpu::spblas::descr_vector_dense & vec_1_descr = data.descr_wvec;
        gpu::spblas::descr_vector_dense & vec_2_descr = (is_factor1_sparse ? data.descr_zvec : data.descr_wvec);
        gpu::spblas::descr_vector_dense & vec_3_descr = ((is_factor1_sparse == is_factor2_sparse) ? data.descr_wvec : data.descr_zvec);

        tm_spmv1.start();
        gpu::spblas::mv<T,I>(hs, 'T', data.descr_Bperm_sp, data.descr_xvec, data.descr_wvec, data.descr_sparse_spmv1, data.buffersize_spmv1, data.buffer_spmv1, 'C');
        if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        tm_spmv1.stop();

        tm_trsv1.start();
        if(do_trsv1_sp_L)  gpu::spblas::trsv<T,I>(hs, 'N', data.descr_L_sp,  vec_1_descr, vec_2_descr, data.descr_sparse_trsv1, data.buffersizes_sptrsv1, data.buffer_sptrsv1_persistent, buffer_tmp, 'C');
        if(do_trsv1_sp_LH) gpu::spblas::trsv<T,I>(hs, 'H', data.descr_LH_sp, vec_1_descr, vec_2_descr, data.descr_sparse_trsv1, data.buffersizes_sptrsv1, data.buffer_sptrsv1_persistent, buffer_tmp, 'C');
        if(do_trsv1_dn_L)  gpu::dnblas::trsv<T,I>(hd, data.n_dofs_domain, data.d_L_dn->vals,  data.ld_domain, 'R', 'N', 'L', vec_1.vals);
        if(do_trsv1_dn_LH) gpu::dnblas::trsv<T,I>(hd, data.n_dofs_domain, data.d_LH_dn->vals, data.ld_domain, 'R', 'H', 'U', vec_1.vals);
        if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        tm_trsv1.stop();

        tm_trsv2.start();
        if(do_trsv2_sp_U)  gpu::spblas::trsv<T,I>(hs, 'N', data.descr_U_sp,  vec_2_descr, vec_3_descr, data.descr_sparse_trsv2, data.buffersizes_sptrsv2, data.buffer_sptrsv2_persistent, buffer_tmp, 'C');
        if(do_trsv2_sp_UH) gpu::spblas::trsv<T,I>(hs, 'H', data.descr_UH_sp, vec_2_descr, vec_3_descr, data.descr_sparse_trsv2, data.buffersizes_sptrsv2, data.buffer_sptrsv2_persistent, buffer_tmp, 'C');
        if(do_trsv2_dn_U)  gpu::dnblas::trsv<T,I>(hd, data.n_dofs_domain, data.d_U_dn->vals,  data.ld_domain, 'R', 'N', 'U', vec_2.vals);
        if(do_trsv2_dn_UH) gpu::dnblas::trsv<T,I>(hd, data.n_dofs_domain, data.d_UH_dn->vals, data.ld_domain, 'R', 'H', 'L', vec_2.vals);
        if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        tm_trsv2.stop();

        tm_spmv2.start();
        gpu::spblas::mv<T,I>(hs, 'N', data.descr_Bperm_sp, vec_3_descr, data.descr_yvec, data.descr_sparse_spmv2, data.buffersize_spmv2, data.buffer_spmv2, 'C');
        if(config->concurrency_apply == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
        tm_spmv2.stop();
    }
    tm_compute.stop();
    
    gpu::dnblas::buffer_unset(hd);

    tm_freeinpool.start();
    {
        gpu::mgm::submit_host_function(q, [&,di,buffer_tmp](){
            domain_data[di].d_L_dn.reset();
            domain_data[di].d_LH_dn.reset();
            domain_data[di].d_U_dn.reset();
            domain_data[di].d_UH_dn.reset();

            void * buffer_tmp_ = buffer_tmp;
            cbmba_res_device->deallocate(buffer_tmp_);
        });
        if(config->concurrency_update == CONCURRENCY::SEQ_WAIT) gpu::mgm::queue_wait(q);
    }
    tm_freeinpool.stop();
}

template <typename T, typename I>
void TotalFETIGpu<T,I>::_apply(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    if(stage != 3) eslog::error("invalid stage when calling apply\n");

    double start = eslog::time();

    gpu::mgm::set_device(device);

    if(is_explicit && config->apply_scatter_gather_where == DEVICE::CPU) apply_explicit_sgcpu(x_cluster, y_cluster);
    if(is_explicit && config->apply_scatter_gather_where == DEVICE::GPU) apply_explicit_sggpu(x_cluster, y_cluster);
    if(is_implicit && config->apply_scatter_gather_where == DEVICE::CPU) apply_implicit_sgcpu(x_cluster, y_cluster);
    if(is_implicit && config->apply_scatter_gather_where == DEVICE::GPU) apply_implicit_sggpu(x_cluster, y_cluster);

    double stop = eslog::time();
    eslog::info("TMP DUAL OPERATOR APPLY TIME:  %12.6f ms\n", (stop - start) * 1000.0);
}

template <typename T, typename I>
void TotalFETIGpu<T,I>::apply(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    apply(x_cluster, y_cluster);
    y_cluster.synchronize();
}

template <typename T, typename I>
void TotalFETIGpu<T,I>::apply(const Matrix_Dual<T> &x_cluster, Matrix_Dual<T> &y_cluster)
{
    Vector_Dual<T> _x, _y;
    _x.size = _y.size = x_cluster.ncols;
    for (int r = 0; r < x_cluster.nrows; ++r) {
        _x.vals = x_cluster.vals + x_cluster.ncols * r;
        _y.vals = y_cluster.vals + y_cluster.ncols * r;
        _apply(_x, _y);
    }
    y_cluster.synchronize();
}

template <typename T, typename I>
void TotalFETIGpu<T,I>::clear_gpu_cache()
{
    size_t size = size_t{1} << 30;
    bool alloc_from_pool = true;
    if(size > cbmba_res_device->get_max_capacity() / 2) alloc_from_pool = false;
    void * data;
    if(alloc_from_pool) {
        cbmba_res_device->do_transaction([&](){ data = cbmba_res_device->allocate(size, align_B); });
    }
    else {
        size_t allocated = 0;
        gpu::mgm::memalloc_device_max(data, allocated, size);
        if(allocated < size) size = allocated;
    }
    gpu::mgm::memset_submit(main_q, data, size, 0);
    gpu::mgm::device_wait();
    if(alloc_from_pool) cbmba_res_device->deallocate(data);
    else gpu::mgm::memfree_device(data);
}

template <typename T, typename I>
void TotalFETIGpu<T,I>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    // just do it on cpu
    #pragma omp parallel for schedule(static,1)
    for (size_t di = 0; di < n_domains; ++di) {
        Vector_Dense<T,I> z;
        z.resize(y[di]);
        applyBt(feti, di, x, z, T{-1});
        math::add(z, T{1}, feti.f[di]);
        domain_data[di].solver_Kreg.solve(z, y[di]);
    }
}

template <typename T, typename I>
void TotalFETIGpu<T,I>::print(const step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/dualop/{Kplus, F}\n"); // Kplus is actually Kreg
        for (size_t di = 0; di < feti.K.size(); ++di) {
            math::store(domain_data[di].Kreg, utils::filename(utils::debugDirectory(step) + "/feti/dualop", (std::string("Kplus") + std::to_string(di)).c_str()).c_str());

            per_domain_stuff & data = domain_data[di];
            Matrix_Dense<T,I,Ah> h_F;
            h_F.resize(data.d_F.nrows, data.d_F.ncols);
            gpu::mgm::copy_submit(main_q, h_F, data.d_F);
            gpu::mgm::queue_wait(main_q);
            if(order_F == 'C') {
                // transpose the matrix to transition colmajor->rowmajor
                for(I r = 0; r < h_F.nrows; r++) {
                    for(I c = 0; c < r; c++) {
                        std::swap(h_F.vals[r * h_F.get_ld() + c], h_F.vals[c * h_F.get_ld() + r]);
                    }
                }
            }
            if(is_system_hermitian) {
                // fill the other triangle
                for(I r = 0; r < h_F.nrows; r++) {
                    I begin = (data.hermitian_F_fill == 'L') ? r+1       : 0;
                    I end   = (data.hermitian_F_fill == 'L') ? h_F.ncols : r;
                    for(I c = begin; c < end; c++) h_F.vals[r * h_F.get_ld() + c] = h_F.vals[c * h_F.get_ld() + r];
                }
            }
            Matrix_Dense<T> F;
            F.shallowCopy(h_F);
            math::store(F, utils::filename(utils::debugDirectory(step) + "/feti/dualop", (std::string("F") + std::to_string(di)).c_str()).c_str());
        }
        math::store(d, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "d").c_str());
    }
}



template<typename U>
void replace_if_auto(U & val, U replace_with)
{
    if(val == U::AUTO) val = replace_with;
}



static TRS2_SOLVE_TYPE use_same_factor(TRS1_SOLVE_TYPE trs1type)
{
    switch(trs1type) {
        case TRS1_SOLVE_TYPE::L: return TRS2_SOLVE_TYPE::UHH;
        case TRS1_SOLVE_TYPE::LHH: return TRS2_SOLVE_TYPE::U;
        default: eslog::error("Invalid trs1type in use_same_factor()\n");
    }
}



template <typename T, typename I>
void TotalFETIGpu<T,I>::config_replace_defaults()
{
    int dimension = info::mesh->dimension;
    [[maybe_unused]] size_t avg_ndofs_per_domain = std::accumulate(feti.K.begin(), feti.K.end(), size_t{0}, [](size_t s, const Matrix_CSR<T,I> & k){ return s + k.nrows; }) / feti.K.size();
    [[maybe_unused]] size_t max_ndofs_per_domain = std::max_element(feti.K.begin(), feti.K.end(), [](const auto & l, const auto & r){ return l.nrows < r.nrows; })->nrows;
    bool is_dense_possible = (max_ndofs_per_domain * max_ndofs_per_domain < (size_t)std::numeric_limits<I>::max());

    if(DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::HERMITIAN_LOWER) {
        std::runtime_error("AUTO gpu options not implemented for lower triangular matrices, todo\n");
        // I could just remove this error (uncomment the if) and use the same options as in the upper triangular case
        //   this would bring a possibly unnecessary transposition, which, if done on the CPU, is hidden behind GPU computations anyway
    }
    if(DirectSparseSolver<T,I>::factorsSymmetry() == Solver_Factors::NONSYMMETRIC_BOTH) {
        std::runtime_error("AUTO gpu options not implemented for unsymmetric systems, todo\n");
    }

    // order of popullating the AUTO options (later can have dependencies on the previous):
    // explicit:
    //   general options:
    //     - F_SHARING_IF_HERMITIAN
    //     - QUEUE_COUNT
    //   set options:
    //     - CONCURRENCY_SET
    //   update options:
    //     - CONCURRENCY_UPDATE
    //     - TRANSPOSE_WHERE
    //     - PATH_IF_HERMITIAN
    //     - TRS1_FACTOR_STORAGE
    //     - TRS2_FACTOR_STORAGE
    //     - TRSM_RHS_SOL_ORDER
    //     - TRS1_SOLVE_TYPE
    //     - TRS2_SOLVE_TYPE
    //   apply options:
    //     - CONCURRENCY_APPLY
    //     - APPLY_SCATTER_GATHER_WHERE
    // implicit:
    //   general+unused options:
    //     - F_SHARING_IF_HERMITIAN
    //     - QUEUE_COUNT
    //     - PATH_IF_HERMITIAN
    //     - TRSM_RHS_SOL_ORDER
    //   set options:
    //     - CONCURRENCY_SET
    //   update options:
    //     - CONCURRENCY_UPDATE
    //   apply options:
    //     - CONCURRENCY_APPLY
    //     - APPLY_SCATTER_GATHER_WHERE
    //     - TRANSPOSE_WHERE
    //     - TRS1_FACTOR_STORAGE
    //     - TRS2_FACTOR_STORAGE
    //     - TRS1_SOLVE_TYPE
    //     - TRS2_SOLVE_TYPE

    if(is_explicit) {
        if(gpu::mgm::get_implementation() == gpu::mgm::gpu_wrapper_impl::CUDA)
        {
            if(gpu::spblas::get_implementation() == gpu::spblas::spblas_wrapper_impl::CUSPARSE_LEGACY) {
                // general
                {
                    replace_if_auto(config->f_sharing_if_hermitian, TRIANGLE_MATRIX_SHARING::SHARED);
                    replace_if_auto(config->queue_count, QUEUE_COUNT::PER_THREAD);
                }

                // set
                {
                    replace_if_auto(config->concurrency_set, CONCURRENCY::PARALLEL);
                }

                // update
                {
                    replace_if_auto(config->concurrency_update, CONCURRENCY::PARALLEL);
                    replace_if_auto(config->transpose_where, DEVICE::CPU);
                    replace_if_auto(config->path_if_hermitian, PATH_IF_HERMITIAN::HERK);

                    MATRIX_STORAGE select_trsm_storage = MATRIX_STORAGE::AUTO;
                    if(dimension == 2)                                  select_trsm_storage = MATRIX_STORAGE::SPARSE;
                    if(dimension == 3 && avg_ndofs_per_domain <  12000) select_trsm_storage = MATRIX_STORAGE::DENSE;
                    if(dimension == 3 && avg_ndofs_per_domain >= 12000) select_trsm_storage = MATRIX_STORAGE::SPARSE;
                    if(!is_dense_possible) select_trsm_storage = MATRIX_STORAGE::SPARSE;
                    replace_if_auto(config->trs1_factor_storage, select_trsm_storage);
                    replace_if_auto(config->trs2_factor_storage, select_trsm_storage);

                    replace_if_auto(config->trsm_rhs_sol_order, MATRIX_ORDER::ROW_MAJOR);

                    TRS1_SOLVE_TYPE select_trsm1_solve_type;
                    bool is_path_herk_used = (config->path_if_hermitian == PATH_IF_HERMITIAN::HERK && DirectSparseSolver<T,I>::factorsSymmetry() != Solver_Factors::NONSYMMETRIC_BOTH);
                    bool is_path_trsm_used = !is_path_herk_used;
                    if(is_path_herk_used && config->trs1_factor_storage == MATRIX_STORAGE::SPARSE) select_trsm1_solve_type = TRS1_SOLVE_TYPE::L;
                    if(is_path_herk_used && config->trs1_factor_storage == MATRIX_STORAGE::DENSE) select_trsm1_solve_type = TRS1_SOLVE_TYPE::LHH;
                    if(is_path_trsm_used) select_trsm1_solve_type = TRS1_SOLVE_TYPE::LHH;
                    replace_if_auto(config->trs1_solve_type, select_trsm1_solve_type);

                    TRS2_SOLVE_TYPE select_trsm2_solve_type = use_same_factor(config->trs1_solve_type);
                    replace_if_auto(config->trs2_solve_type, select_trsm2_solve_type);
                }

                // apply
                {
                    replace_if_auto(config->concurrency_apply, CONCURRENCY::PARALLEL);
                    replace_if_auto(config->apply_scatter_gather_where, DEVICE::GPU);
                }
            }
            else if(gpu::spblas::get_implementation() == gpu::spblas::spblas_wrapper_impl::CUSPARSE_MODERN) {
                // general
                {
                    replace_if_auto(config->f_sharing_if_hermitian, TRIANGLE_MATRIX_SHARING::SHARED);
                    replace_if_auto(config->queue_count, QUEUE_COUNT::PER_THREAD);
                }

                // set
                {
                    replace_if_auto(config->concurrency_set, CONCURRENCY::PARALLEL);
                }

                // update
                {
                    replace_if_auto(config->concurrency_update, CONCURRENCY::PARALLEL);
                    replace_if_auto(config->transpose_where, DEVICE::CPU);
                    replace_if_auto(config->path_if_hermitian, PATH_IF_HERMITIAN::HERK);

                    MATRIX_STORAGE select_trsm_storage = MATRIX_STORAGE::DENSE;
                    if(!is_dense_possible) select_trsm_storage = MATRIX_STORAGE::SPARSE;
                    replace_if_auto(config->trs1_factor_storage, select_trsm_storage);
                    replace_if_auto(config->trs2_factor_storage, select_trsm_storage);

                    MATRIX_ORDER select_order = MATRIX_ORDER::AUTO;
                    if(dimension == 2) select_order = MATRIX_ORDER::COL_MAJOR;
                    if(dimension == 3) select_order = MATRIX_ORDER::ROW_MAJOR;
                    replace_if_auto(config->trsm_rhs_sol_order, select_order);

                    replace_if_auto(config->trs1_solve_type, TRS1_SOLVE_TYPE::LHH);
                    replace_if_auto(config->trs2_solve_type, TRS2_SOLVE_TYPE::U);
                }

                // apply
                {
                    replace_if_auto(config->concurrency_apply, CONCURRENCY::PARALLEL);
                    replace_if_auto(config->apply_scatter_gather_where, DEVICE::GPU);
                }
            }
            else {
                eslog::error("Unexpected gpu sparse blas implementation\n");
            }
        }
        if(gpu::mgm::get_implementation() == gpu::mgm::gpu_wrapper_impl::ROCM) {
            // general
            {
                replace_if_auto(config->f_sharing_if_hermitian, TRIANGLE_MATRIX_SHARING::SHARED);
                replace_if_auto(config->queue_count, QUEUE_COUNT::PER_THREAD);
            }

            // set
            {
                replace_if_auto(config->concurrency_set, CONCURRENCY::PARALLEL);
            }

            // update
            {
                replace_if_auto(config->concurrency_update, CONCURRENCY::PARALLEL);
                replace_if_auto(config->transpose_where, DEVICE::CPU);
                replace_if_auto(config->path_if_hermitian, PATH_IF_HERMITIAN::HERK);

                MATRIX_STORAGE select_trsm1_storage = MATRIX_STORAGE::SPARSE;
                if(dimension == 3 && avg_ndofs_per_domain < 32000) select_trsm1_storage = MATRIX_STORAGE::DENSE;
                if(!is_dense_possible) select_trsm1_storage = MATRIX_STORAGE::SPARSE;
                replace_if_auto(config->trs1_factor_storage, select_trsm1_storage);

                MATRIX_STORAGE select_trsm2_storage = MATRIX_STORAGE::DENSE;
                if(dimension == 2 && avg_ndofs_per_domain > 46000) select_trsm2_storage = MATRIX_STORAGE::SPARSE;
                if(!is_dense_possible) select_trsm2_storage = MATRIX_STORAGE::SPARSE;
                replace_if_auto(config->trs2_factor_storage, select_trsm2_storage);

                replace_if_auto(config->trsm_rhs_sol_order, MATRIX_ORDER::ROW_MAJOR);

                TRS1_SOLVE_TYPE select_trsm1_solve_type = TRS1_SOLVE_TYPE::L;
                bool is_path_herk_used = (config->path_if_hermitian == PATH_IF_HERMITIAN::HERK && DirectSparseSolver<T,I>::factorsSymmetry() != Solver_Factors::NONSYMMETRIC_BOTH);
                bool is_path_trsm_used = !is_path_herk_used;
                if(is_path_trsm_used && config->trs1_factor_storage == MATRIX_STORAGE::SPARSE) select_trsm1_solve_type = TRS1_SOLVE_TYPE::LHH;
                replace_if_auto(config->trs1_solve_type, select_trsm1_solve_type);

                replace_if_auto(config->trs2_solve_type, TRS2_SOLVE_TYPE::U);
            }

            // apply
            {
                replace_if_auto(config->concurrency_apply, CONCURRENCY::PARALLEL);
                replace_if_auto(config->apply_scatter_gather_where, DEVICE::GPU);
            }
        }
        if(gpu::mgm::get_implementation() == gpu::mgm::gpu_wrapper_impl::ONEAPI) {
            eslog::error("not implemented\n");
        }
    }
    if(is_implicit) {
        if(gpu::mgm::get_implementation() == gpu::mgm::gpu_wrapper_impl::CUDA)
        {
            if(gpu::spblas::get_implementation() == gpu::spblas::spblas_wrapper_impl::CUSPARSE_LEGACY) {
                // general+unused
                {
                    // replace_if_auto(config->f_sharing_if_hermitian, TRIANGLE_MATRIX_SHARING::SHARED);
                    replace_if_auto(config->queue_count, QUEUE_COUNT::PER_THREAD);
                    // replace_if_auto(config->path_if_hermitian, PATH_IF_HERMITIAN::TRSM);
                    // replace_if_auto(config->trsm_rhs_sol_order, MATRIX_ORDER::ROW_MAJOR);
                }

                // set
                {
                    replace_if_auto(config->concurrency_set, CONCURRENCY::PARALLEL);
                }

                // update
                {
                    replace_if_auto(config->concurrency_update, CONCURRENCY::PARALLEL);
                }

                // apply
                {
                    replace_if_auto(config->concurrency_apply, CONCURRENCY::PARALLEL);
                    replace_if_auto(config->apply_scatter_gather_where, DEVICE::GPU);
                    replace_if_auto(config->transpose_where, DEVICE::CPU);

                    MATRIX_STORAGE select_trsv1_storage = MATRIX_STORAGE::DENSE;
                    if(dimension == 2 && avg_ndofs_per_domain > 8000) select_trsv1_storage = MATRIX_STORAGE::SPARSE;
                    if(!is_dense_possible) select_trsv1_storage = MATRIX_STORAGE::SPARSE;
                    replace_if_auto(config->trs1_factor_storage, select_trsv1_storage);

                    MATRIX_STORAGE select_trsv2_storage = MATRIX_STORAGE::DENSE;
                    if(dimension == 2 && avg_ndofs_per_domain > 23000) select_trsv2_storage = MATRIX_STORAGE::SPARSE;
                    if(!is_dense_possible) select_trsv2_storage = MATRIX_STORAGE::SPARSE;
                    replace_if_auto(config->trs2_factor_storage, select_trsv2_storage);

                    replace_if_auto(config->trs1_solve_type, TRS1_SOLVE_TYPE::L);

                    TRS2_SOLVE_TYPE select_trsv2_solve_type = TRS2_SOLVE_TYPE::AUTO;
                    if(config->trs2_factor_storage == MATRIX_STORAGE::SPARSE) select_trsv2_solve_type = TRS2_SOLVE_TYPE::U;
                    if(config->trs2_factor_storage == MATRIX_STORAGE::DENSE)  select_trsv2_solve_type = TRS2_SOLVE_TYPE::UHH;
                    replace_if_auto(config->trs2_solve_type, select_trsv2_solve_type);
                }
            }
            else if(gpu::spblas::get_implementation() == gpu::spblas::spblas_wrapper_impl::CUSPARSE_MODERN) {
                // general+unused
                {
                    // replace_if_auto(config->f_sharing_if_hermitian, TRIANGLE_MATRIX_SHARING::SHARED);
                    replace_if_auto(config->queue_count, QUEUE_COUNT::PER_THREAD);
                    // replace_if_auto(config->path_if_hermitian, PATH_IF_HERMITIAN::TRSM);
                    // replace_if_auto(config->trsm_rhs_sol_order, MATRIX_ORDER::ROW_MAJOR);
                }

                // set
                {
                    replace_if_auto(config->concurrency_set, CONCURRENCY::PARALLEL);
                }

                // update
                {
                    replace_if_auto(config->concurrency_update, CONCURRENCY::PARALLEL);
                }

                // apply
                {
                    replace_if_auto(config->concurrency_apply, CONCURRENCY::SEQ_CONTINUE);
                    replace_if_auto(config->apply_scatter_gather_where, DEVICE::GPU);
                    replace_if_auto(config->transpose_where, DEVICE::CPU);

                    MATRIX_STORAGE select_trsv_storage = MATRIX_STORAGE::SPARSE;
                    if(dimension == 3 && avg_ndofs_per_domain < 8000) select_trsv_storage = MATRIX_STORAGE::DENSE;
                    if(!is_dense_possible) select_trsv_storage = MATRIX_STORAGE::SPARSE;
                    replace_if_auto(config->trs1_factor_storage, select_trsv_storage);
                    replace_if_auto(config->trs2_factor_storage, select_trsv_storage);

                    replace_if_auto(config->trs1_solve_type, TRS1_SOLVE_TYPE::L);
                    replace_if_auto(config->trs2_solve_type, TRS2_SOLVE_TYPE::UHH);
                }
            }
            else {
                eslog::error("Unexpected gpu sparse blas implementation\n");
            }
        }
        if(gpu::mgm::get_implementation() == gpu::mgm::gpu_wrapper_impl::ROCM) {
            // general+unused
            {
                // replace_if_auto(config->f_sharing_if_hermitian, TRIANGLE_MATRIX_SHARING::SHARED);
                replace_if_auto(config->queue_count, QUEUE_COUNT::PER_THREAD);
                // replace_if_auto(config->path_if_hermitian, PATH_IF_HERMITIAN::TRSM);
                // replace_if_auto(config->trsm_rhs_sol_order, MATRIX_ORDER::ROW_MAJOR);
            }

            // set
            {
                replace_if_auto(config->concurrency_set, CONCURRENCY::PARALLEL);
            }

            // update
            {
                replace_if_auto(config->concurrency_update, CONCURRENCY::PARALLEL);
            }

            // apply
            {
                replace_if_auto(config->concurrency_apply, CONCURRENCY::SEQ_CONTINUE);
                replace_if_auto(config->apply_scatter_gather_where, DEVICE::GPU);
                replace_if_auto(config->transpose_where, DEVICE::CPU);

                MATRIX_STORAGE select_trsv_storage = MATRIX_STORAGE::SPARSE;
                if(dimension == 3 && avg_ndofs_per_domain < 16000) select_trsv_storage = MATRIX_STORAGE::DENSE;
                if(!is_dense_possible) select_trsv_storage = MATRIX_STORAGE::SPARSE;
                replace_if_auto(config->trs1_factor_storage, select_trsv_storage);
                replace_if_auto(config->trs2_factor_storage, select_trsv_storage);

                TRS1_SOLVE_TYPE select_trsv1_solve_type = TRS1_SOLVE_TYPE::AUTO;
                if(config->trs1_factor_storage == MATRIX_STORAGE::SPARSE) select_trsv1_solve_type = TRS1_SOLVE_TYPE::L;
                if(config->trs1_factor_storage == MATRIX_STORAGE::DENSE)  select_trsv1_solve_type = TRS1_SOLVE_TYPE::LHH;
                replace_if_auto(config->trs1_solve_type, select_trsv1_solve_type);

                replace_if_auto(config->trs2_solve_type, TRS2_SOLVE_TYPE::U);
            }
        }
        if(gpu::mgm::get_implementation() == gpu::mgm::gpu_wrapper_impl::ONEAPI) {
            eslog::error("not implemented\n");
        }
    }
}



#define INSTANTIATE_T_I(T,I) \
template class TotalFETIGpu<T,I>;

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
