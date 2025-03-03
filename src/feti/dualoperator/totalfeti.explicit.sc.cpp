
#include "totalfeti.explicit.sc.h"
#include "math/wrappers/math.blas.h"
#include "feti/common/applyB.h"
#include "basis/utilities/minmaxavg.h"
#include "my_timer.h"
#include "gpu/gpu_kernels.h"

#include <algorithm>



namespace espreso {

template<typename T, typename I>
TotalFETIExplicitSc<T,I>::TotalFETIExplicitSc(FETI<T> &feti, bool apply_on_gpu_)
: DualOperator<T>(feti), apply_on_gpu(apply_on_gpu_)
{
    if(apply_on_gpu && !gpu::mgm::is_linked()) {
        eslog::error("TotalFETIExplicitSc: apply on GPU not supported, espreso compiled without GPU support\n");
    }
}



template<typename T, typename I>
TotalFETIExplicitSc<T,I>::~TotalFETIExplicitSc()
{
    // my_timer tm_total;

    // tm_total.start();
    if(apply_on_gpu) {
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_stuff & data = domain_data[di];
            data.d_F.clear();
            data.d_applyg_D2C.clear();
            data.d_apply_x.clear();
            data.d_apply_y.clear();
            data.d_apply_z.clear();
            data.d_apply_w.clear();
        }

        for(size_t qi = 0; qi < n_queues; qi++) {
            gpu::dnblas::buffer_unset(handles_dense[qi]);
            gpu::mgm::memfree_device(d_buffers_dense[qi]);
        }

        d_Fs_allocated.clear();
        d_applyg_x_cluster.clear();
        d_applyg_y_cluster.clear();
        d_applyg_xs_pointers.clear();
        d_applyg_ys_pointers.clear();
        d_applyg_n_dofs_interfaces.clear();
        d_applyg_D2Cs_pointers.clear();

        gpu::mgm::queue_destroy(main_q);
        for(gpu::mgm::queue & q : queues) gpu::mgm::queue_destroy(q);
        for(gpu::dnblas::handle & h : handles_dense) gpu::dnblas::handle_destroy(h);

        d_Fs_allocated.clear();
    }
    // tm_total.stop();

    // print_timer("Destroy total", tm_total);
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::info()
{
    eslog::info(" = EXPLICIT TOTAL FETI OPERATOR USING SPBLAS SCHUR COMPLEMENT                                = \n");
    eslog::info(" =   EXTERNAL SCHUL COMPLEMENT SOLVER     %50s = \n", SchurComplementSolver<T,I>::name());
    eslog::info(" =   %-50s       %+30s = \n", "APPLY_WHERE", apply_on_gpu ? "GPU" : "CPU");
    eslog::info(minmaxavg<double>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.F.nrows * data.F.get_ld() * sizeof(T) / (1024.0 * 1024.0); }).to_string("  F MEMORY [MB]").c_str());
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_domain; }).to_string("  Domain volume [dofs]").c_str());
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_stuff & data){ return data.n_dofs_interface; }).to_string("  Domain surface [dofs]").c_str());
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::set(const step::Step &step)
{
    gpu_wait_update = false;
    gpu_wait_apply = false;
    gpu_apply_parallel = (gpu::mgm::get_implementation() != gpu::mgm::gpu_wrapper_impl::ONEAPI);
    is_system_hermitian = std::all_of(feti.K.begin(), feti.K.end(), [](Matrix_CSR<T,I> & K){ return getSymmetry(K.type) == Matrix_Symmetry::HERMITIAN; });
    Fs_share_memory = is_system_hermitian;



    // my_timer tm_total, tm_gpuinit, tm_gpuset, tm_gpucreate, tm_alloc_Fs, tm_mainloop, tm_inner, tm_prepare, tm_commit, tm_symbfact, tm_gpualloc, tm_gpuapplystuff;

    if(apply_on_gpu) {
        device = gpu::mgm::get_device_by_mpi(info::mpi::rank, info::mpi::size);
    }

    // tm_total.start();

    n_domains = feti.K.size();
    domain_data.resize(n_domains);
    for(size_t di = 0; di < n_domains; di++) {
        domain_data[di].n_dofs_interface = feti.B1[di].nrows;
        domain_data[di].n_dofs_domain = feti.B1[di].ncols;
    }

    if(apply_on_gpu) {
        // tm_gpuinit.start();
        gpu::mgm::init_gpu(device);
        // tm_gpuinit.stop();

        // tm_gpuset.start();
        gpu::mgm::set_device(device);
        // tm_gpuset.stop();

        // tm_gpucreate.start();
        n_queues = omp_get_max_threads();
        queues.resize(n_queues);
        handles_dense.resize(n_queues);
        d_buffers_dense.resize(n_queues);
        gpu::mgm::queue_create(main_q);
        for(gpu::mgm::queue & q : queues) gpu::mgm::queue_create(q);
        for(size_t i = 0; i < n_queues; i++) gpu::dnblas::handle_create(handles_dense[i], queues[i]);
        // tm_gpucreate.stop();
    }

    // tm_alloc_Fs.start();
    if(Fs_share_memory) {
        Fs_allocated.resize((n_domains - 1) / 2 + 1);
        std::vector<size_t> domain_idxs_sorted_by_f_size_desc(n_domains);
        for(size_t di = 0; di < n_domains; di++) {
            domain_idxs_sorted_by_f_size_desc[di] = di;
        }
        std::sort(domain_idxs_sorted_by_f_size_desc.rbegin(), domain_idxs_sorted_by_f_size_desc.rend(), [&](size_t dl, size_t dr){ return domain_data[dl].n_dofs_interface < domain_data[dr].n_dofs_interface; });
        for(size_t i = 0; i < n_domains; i++) {
            size_t di = domain_idxs_sorted_by_f_size_desc[i];
            auto & data = domain_data[di];

            if(i % 2 == 0) {
                Fs_allocated[i / 2].resize(data.n_dofs_interface + 1, data.n_dofs_interface);
                data.F_fill = 'U';
                data.F.shallowCopy(Fs_allocated[i / 2]);
                data.F.nrows = data.n_dofs_interface;
                data.F.ncols = data.n_dofs_interface;
            }
            else {
                data.F_fill = 'L';
                data.F.shallowCopy(Fs_allocated[i / 2]);
                data.F.nrows = data.n_dofs_interface;
                data.F.ncols = data.n_dofs_interface;
                data.F.vals += data.F.get_ld();
            }
            if constexpr(utils::is_real<T>())    data.F.type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE;
            if constexpr(utils::is_complex<T>()) data.F.type = Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE;
        }

        if(apply_on_gpu) {
            d_Fs_allocated.resize((n_domains - 1) / 2 + 1);
            for(size_t i = 0; i < n_domains; i++) {
                size_t di = domain_idxs_sorted_by_f_size_desc[i];
                auto & data = domain_data[di];

                if(i % 2 == 0) {
                    d_Fs_allocated[i / 2].resize(data.n_dofs_interface + 1, data.n_dofs_interface);
                    // data.F_fill = 'U';
                    data.d_F.shallowCopy(d_Fs_allocated[i / 2]);
                    data.d_F.nrows = data.n_dofs_interface;
                    data.d_F.ncols = data.n_dofs_interface;
                }
                else {
                    // data.F_fill = 'L';
                    data.d_F.shallowCopy(d_Fs_allocated[i / 2]);
                    data.d_F.nrows = data.n_dofs_interface;
                    data.d_F.ncols = data.n_dofs_interface;
                    data.d_F.vals += data.d_F.get_ld();
                }
                if constexpr(utils::is_real<T>())    data.d_F.type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE;
                if constexpr(utils::is_complex<T>()) data.d_F.type = Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE;
            }
        }
    }
    else {
        for(size_t di = 0; di < n_domains; di++) {
            auto & data = domain_data[di];

            data.F_fill = 'U';
            data.F.resize(data.n_dofs_interface, data.n_dofs_interface, data.n_dofs_interface);
            if constexpr(utils::is_real<T>())    data.F.type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE;
            if constexpr(utils::is_complex<T>()) data.F.type = Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE;
        }

        if(apply_on_gpu) {
            for(size_t di = 0; di < n_domains; di++) {
                auto & data = domain_data[di];

                // data.F_fill = 'U';
                data.d_F.resize(data.n_dofs_interface, data.n_dofs_interface, data.n_dofs_interface);
                if constexpr(utils::is_real<T>())    data.d_F.type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE;
                if constexpr(utils::is_complex<T>()) data.d_F.type = Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE;
            }
        }
    }
    // tm_alloc_Fs.stop();

    // tm_mainloop.start();
    #pragma omp parallel for schedule(static,1)
    for(size_t di = 0; di < n_domains; di++) {
        // tm_inner.start();
        auto & data = domain_data[di];

        if(getSymmetry(feti.K[di].type) != Matrix_Symmetry::HERMITIAN || feti.K[di].shape != Matrix_Shape::UPPER) {
            eslog::error("implemented only for hermitian K matrices stored in upper triangle. TODO\n");
        }

        // tm_prepare.start();
        data.Kreg.type = feti.K[di].type;
        data.Kreg.shape = feti.K[di].shape;
        math::combine(data.Kreg, feti.K[di], feti.RegMat[di]);

        Matrix_CSR<T,I> & B = feti.B1[di];
        data.Bt.resize(B.ncols, B.nrows, B.nnz);
        if constexpr(utils::is_real<T>())    data.Bt.type = Matrix_Type::REAL_NONSYMMETRIC;
        if constexpr(utils::is_complex<T>()) data.Bt.type = Matrix_Type::COMPLEX_NONSYMMETRIC;
        data.map_B_transpose.resize(B.nnz);
        math::csrTranspose(data.Bt, B, data.map_B_transpose, math::CsrTransposeStage::Pattern, true);

        data.null_matrix_A21.resize(data.n_dofs_interface, data.n_dofs_domain, 0);
        std::fill_n(data.null_matrix_A21.rows, data.null_matrix_A21.nrows + 1, 0);
        if constexpr(utils::is_real<T>())    data.null_matrix_A21.type = Matrix_Type::REAL_NONSYMMETRIC;
        if constexpr(utils::is_complex<T>()) data.null_matrix_A21.type = Matrix_Type::COMPLEX_NONSYMMETRIC;
        data.null_matrix_A22.resize(data.n_dofs_interface, data.n_dofs_interface, 0);
        std::fill_n(data.null_matrix_A22.rows, data.null_matrix_A22.nrows + 1, 0);
        if constexpr(utils::is_real<T>())    data.null_matrix_A22.type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE;
        if constexpr(utils::is_complex<T>()) data.null_matrix_A22.type = Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE;
        data.null_matrix_A22.shape = Matrix_Shape::UPPER;
        // tm_prepare.stop();

        // tm_commit.start();
        data.sc_solver.commitMatrix(data.Kreg, data.Bt, data.null_matrix_A21, data.null_matrix_A22);
        // tm_commit.stop();
        // tm_symbfact.start();
        data.sc_solver.factorizeSymbolic();
        // tm_symbfact.stop();

        if(apply_on_gpu) {
            gpu::dnblas::handle & hd = handles_dense[di % n_queues];
            gpu::dnblas::buffer_collect_size(hd, data.buffersize_tmp_apply, [&](){
                T * dummyptrT = reinterpret_cast<T*>(sizeof(T));
                if( is_system_hermitian) gpu::dnblas::hemv(hd, data.n_dofs_interface, dummyptrT, data.d_F.get_ld(), 'R', 'N', data.F_fill, dummyptrT, dummyptrT);
                if(!is_system_hermitian) gpu::dnblas::gemv(hd, data.n_dofs_interface, data.n_dofs_interface, dummyptrT, data.d_F.get_ld(), 'R', 'N', dummyptrT, dummyptrT);
            });

            // tm_gpualloc.start();
            data.d_applyg_D2C.resize(data.n_dofs_interface);
            data.d_apply_x.resize(data.n_dofs_interface);
            data.d_apply_y.resize(data.n_dofs_interface);
            data.d_apply_z.resize(data.n_dofs_domain);
            data.d_apply_w.resize(data.n_dofs_domain);
            // tm_gpualloc.stop();
        }
        else {
            data.x.resize(data.n_dofs_interface);
            data.y.resize(data.n_dofs_interface);
        }
        // tm_inner.stop();
    }
    // tm_mainloop.stop();

    if(apply_on_gpu) {
        for(size_t qi = 0; qi < n_queues; qi++) {
            size_t buffersize_tmp_apply_max = 0;
            for(size_t di = qi; di < n_domains; di += n_queues) {
                buffersize_tmp_apply_max = std::max(buffersize_tmp_apply_max, domain_data[di].buffersize_tmp_apply);
            }
            d_buffers_dense[qi] = gpu::mgm::memalloc_device(buffersize_tmp_apply_max);
            gpu::dnblas::buffer_set(handles_dense[qi], d_buffers_dense[qi], buffersize_tmp_apply_max);
        }
        // tm_gpuapplystuff.start();
        d_applyg_x_cluster.resize(feti.lambdas.size);
        d_applyg_y_cluster.resize(feti.lambdas.size);
        Vector_Dense<T*,I,gpu::mgm::Ah> h_applyg_xs_pointers;
        Vector_Dense<T*,I,gpu::mgm::Ah> h_applyg_ys_pointers;
        Vector_Dense<I,I,gpu::mgm::Ah> h_applyg_n_dofs_interfaces;
        Vector_Dense<I*,I,gpu::mgm::Ah> h_applyg_D2Cs_pointers;
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
        gpu::mgm::queue_wait(main_q);
        // tm_gpuapplystuff.stop();
    }

    // tm_total.stop();

    // print_timer("Set     total", tm_total);
    // print_timer("Set       gpuinit", tm_gpuinit);
    // print_timer("Set       gpuset", tm_gpuset);
    // print_timer("Set       gpucreate", tm_gpucreate);
    // print_timer("Set       alloc_Fs", tm_alloc_Fs);
    // print_timer("Set       mainloop", tm_mainloop);
    // print_timer("Set         inner", tm_inner);
    // print_timer("Set           prepare", tm_prepare);
    // print_timer("Set           commit", tm_commit);
    // print_timer("Set           symbfact", tm_symbfact);
    // print_timer("Set           gpualloc", tm_gpualloc);
    // print_timer("Set       gpuapplystuff", tm_gpuapplystuff);
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::update(const step::Step &step)
{
    // my_timer tm_total, tm_assemble, tm_inner, tm_updatevals, tm_factnumer, tm_copyF, tm_maked, tm_wait;

    // tm_total.start();
    // tm_assemble.start();
    #pragma omp parallel for schedule(static,1)
    for(size_t di = 0; di < n_domains; di++) {
        // tm_inner.start();
        auto & data = domain_data[di];

        math::sumCombined(data.Kreg, T{1.0}, feti.K[di], feti.RegMat[di]);

        Matrix_CSR<T,I> & B = feti.B1[di];
        math::csrTranspose(data.Bt, B, data.map_B_transpose, math::CsrTransposeStage::Values, true);

        // tm_updatevals.start();
        data.sc_solver.updateMatrixValues();
        // tm_updatevals.stop();
        // tm_factnumer.start();
        data.sc_solver.factorizeNumericAndGetSc(data.F, data.F_fill, T{-1});
        // tm_factnumer.stop();
        // tm_inner.stop();
    }
    // tm_assemble.stop();

    if(apply_on_gpu) {
        // tm_copyF.start();
        for(size_t i = 0; i < d_Fs_allocated.size(); i++) {
            gpu::mgm::copy_submit(main_q, d_Fs_allocated[i], Fs_allocated[i]);
        }
        if(gpu_wait_update) gpu::mgm::queue_wait(main_q);
        // tm_copyF.stop();
    }
    // tm_total.stop();

    // tm_maked.start();
    {
        if (feti.updated.B) {
            d.resize();
        }
        std::vector<Vector_Dense<T,I>> Kplus_fs(n_domains);
        #pragma omp parallel for schedule(static,1)
        for(size_t di = 0; di < n_domains; di++) {
            auto & data = domain_data[di];
            Kplus_fs[di].resize(data.n_dofs_domain);
            data.sc_solver.solveA11(feti.f[di], Kplus_fs[di]);
        }
        applyB(feti, Kplus_fs, d);
        d.synchronize();
        math::add(d, T{-1}, feti.c);
    }
    // tm_maked.stop();

    // tm_wait.start();
    if(apply_on_gpu) {
        gpu::mgm::queue_wait(main_q);
    }
    // tm_wait.stop();

    // print_timer("Update  total", tm_total);
    // print_timer("Update    assemble", tm_assemble);
    // print_timer("Update      inner", tm_inner);
    // print_timer("Update        updatevals", tm_updatevals);
    // print_timer("Update        factnumer", tm_factnumer);
    // print_timer("Update    copyF", tm_copyF);
    // print_timer("Update  maked", tm_maked);
    // print_timer("Update  wait", tm_wait);
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::_apply_cpu(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    // my_timer tm_total, tm_mv;

    // tm_total.start();
    memset(y_cluster.vals, 0, y_cluster.size * sizeof(T));

    // #pragma omp parallel for schedule(static,1)
    for(size_t di = 0; di < n_domains; di++) {
        auto & data = domain_data[di];

        std::vector<I> & D2C = feti.D2C[di];

        for(I i = 0; i < data.n_dofs_interface; i++) {
            data.x.vals[i] = x_cluster.vals[D2C[i]];
        }

        // tm_mv.start();
        if( is_system_hermitian) math::blas::apply_hermitian(data.y, T{1}, data.F, data.F_fill, T{0}, data.x);
        if(!is_system_hermitian) eslog::error("not implemented for non-hermitian\n");
        // tm_mv.stop();

        for(I i = 0; i < data.n_dofs_interface; i++) {
            #pragma omp atomic
            y_cluster.vals[D2C[i]] += data.y.vals[i];
        }
    }
    // tm_total.stop();

    // print_timer("ApplyC  total", tm_total);
    // print_timer("ApplyC    mv", tm_mv);
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::_apply_gpu(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    // my_timer tm_total;
    // my_timer tm_copyin, tm_scatter, tm_mv_outer, tm_mv, tm_zerofill, tm_gather, tm_copyout, tm_wait;

    // tm_total.start();

    // copy x_cluster to device
    // tm_copyin.start();
    gpu::mgm::copy_submit(main_q, d_applyg_x_cluster, x_cluster);
    if(gpu_wait_apply) gpu::mgm::queue_wait(main_q);
    // tm_copyin.stop();

    // scatter
    // tm_scatter.start();
    gpu::kernels::DCmap_scatter(main_q, d_applyg_xs_pointers, d_applyg_n_dofs_interfaces, d_applyg_x_cluster, d_applyg_D2Cs_pointers);
    if(gpu_wait_apply) gpu::mgm::queue_wait(main_q);
    // tm_scatter.stop();

    gpu::mgm::queue_async_barrier({main_q}, queues);

    // apply
    // tm_mv_outer.start();
    #pragma omp parallel for schedule(static,1) if(gpu_apply_parallel)
    for(size_t di = 0; di < n_domains; di++) {
        gpu::mgm::queue & q = queues[di % n_queues];
        gpu::dnblas::handle & hd = handles_dense[di % n_queues];
        per_domain_stuff & data = domain_data[di];

        // tm_mv.start();
        if( is_system_hermitian) gpu::dnblas::hemv<T,I>(hd, data.d_F.nrows, data.d_F.vals, data.d_F.get_ld(), 'R', 'N', data.F_fill, data.d_apply_x.vals, data.d_apply_y.vals);
        if(!is_system_hermitian) gpu::dnblas::gemv<T,I>(hd, data.d_F.nrows, data.d_F.ncols, data.d_F.vals, data.d_F.get_ld(), 'R', 'N', data.d_apply_x.vals, data.d_apply_y.vals);
        if(gpu_wait_apply) gpu::mgm::queue_wait(q);
        // tm_mv.stop();
    }
    // tm_mv_outer.stop();

    // zerofill y_cluster on device
    // tm_zerofill.start();
    gpu::mgm::memset_submit(main_q, d_applyg_y_cluster.vals, d_applyg_y_cluster.size * sizeof(T), 0);
    if(gpu_wait_apply) gpu::mgm::queue_wait(main_q);
    // tm_zerofill.stop();

    gpu::mgm::queue_async_barrier(queues, {main_q});

    // gather
    // tm_gather.start();
    gpu::kernels::DCmap_gather(main_q, d_applyg_ys_pointers, d_applyg_n_dofs_interfaces, d_applyg_y_cluster, d_applyg_D2Cs_pointers);
    if(gpu_wait_apply) gpu::mgm::queue_wait(main_q);
    // tm_gather.stop();

    // copy y_cluster from device
    // tm_copyout.start();
    gpu::mgm::copy_submit(main_q, y_cluster, d_applyg_y_cluster);
    if(gpu_wait_apply) gpu::mgm::queue_wait(main_q);
    // tm_copyout.stop();

    // wait
    // tm_wait.start();
    gpu::mgm::device_wait();
    // tm_wait.stop();

    // tm_total.stop();

    // print_timer("ApplyG  total", tm_total);
    // print_timer("ApplyG    copyin", tm_copyin);
    // print_timer("ApplyG    scatter", tm_scatter);
    // print_timer("ApplyG    mv_outer", tm_mv_outer);
    // print_timer("ApplyG      mv", tm_mv);
    // print_timer("ApplyG    zerofill", tm_zerofill);
    // print_timer("ApplyG    gather", tm_gather);
    // print_timer("ApplyG    copyout", tm_copyout);
    // print_timer("ApplyG    wait", tm_wait);
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::_apply(const Vector_Dual<T> &x_cluster, Vector_Dual<T> &y_cluster)
{
    if(apply_on_gpu) {
        _apply_gpu(x_cluster, y_cluster);
    }
    else {
        _apply_cpu(x_cluster, y_cluster);
    }
}



template <typename T, typename I>
void TotalFETIExplicitSc<T,I>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    _apply(x, y);
    y.synchronize();
}



template <typename T, typename I>
void TotalFETIExplicitSc<T,I>::apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y)
{
    Vector_Dual<T> _x, _y;
    _x.size = _y.size = x.ncols;
    for (int r = 0; r < x.nrows; ++r) {
        _x.vals = x.vals + x.ncols * r;
        _y.vals = y.vals + y.ncols * r;
        _apply(_x, _y);
    }
    y.synchronize();
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    #pragma omp parallel for schedule(static,1)
    for (size_t di = 0; di < n_domains; ++di) {
        Vector_Dense<T,I> z;
        z.resize(y[di]);
        applyBt(feti, di, x, z, T{-1});
        math::add(z, T{1}, feti.f[di]);
        domain_data[di].sc_solver.solveA11(z, y[di]);
    }
}



template<typename T, typename I>
void TotalFETIExplicitSc<T,I>::print(const step::Step &step)
{
    eslog::error("TotalFETIExplicitSc::print not implemented");
}





#define INSTANTIATE_T_I(T,I) \
template class TotalFETIExplicitSc<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        /* INSTANTIATE_T(std::complex<double>) */

    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I

}
