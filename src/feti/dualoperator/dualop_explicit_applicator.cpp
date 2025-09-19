
#include "dualop_explicit_applicator.h"
#include "math/primitives_new/allocator_new.h"
#include "gpu/gpu_kernels.h"
#include "math/wrappers/math.blas.h"
#include "math/operations/fill_dnx.h"
#include "math/operations/submatrix_dnx_dnx_noncontig.h"
#include "math/operations/supermatrix_dnx_dnx_noncontig.h"
#include "gpu/operations/submatrix_ddnx_ddnx_noncontig.h"
#include "gpu/operations/supermatrix_ddnx_ddnx_noncontig.h"
#include "gpu/operations/hemm_ddnx_ddny_ddnz.h"
#include "gpu/operations/gemm_ddnx_ddny_ddnz.h"
#include "basis/utilities/stacktimer.h"

#include <numeric>



namespace espreso {



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::set_config(bool wait_intermediate_, bool timers_inner_)
{
    wait_intermediate = wait_intermediate_;
    timers_inner = timers_inner_;
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::set_handles(gpu::mgm::queue * main_q_, std::vector<gpu::mgm::queue> * queues_, std::vector<gpu::dnblas::handle> * handles_dnblas_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    main_q = main_q_;
    queues = queues_;
    handles_dnblas = handles_dnblas_;

    if(main_q == nullptr || queues == nullptr || handles_dnblas == nullptr) eslog::error("wrong handles\n");

    called_set_handles = true;
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::set_feti(FETI<T> * feti)
{
    n_dofs_cluster_interface = feti->lambdas.size;
    use_gpu = feti->use_gpu;
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::set_memory(char vector_mem_, char Fs_mem_)
{
    if(called_setup) eslog::error("cannot change vector memory anymore\n");

    vector_mem = vector_mem_;
    Fs_mem = Fs_mem_;
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::set_D2C_map(std::vector<std::vector<I>> * D2C_old_)
{
    if(D2C_old != nullptr) eslog::error("D2C map has already been set\n");

    D2C_old = D2C_old_;
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::set_Fs(std::vector<MatrixDenseView_new<T>*> Fs_)
{
    if(Fs.size() > 0) eslog::error("Fs have already been set\n");

    Fs = std::move(Fs_);
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::set_apply_target(char target_)
{
    if(called_setup) eslog::error("cannot change apply target\n");

    apply_target = target_;
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::setup()
{
    stacktimer::push("dualop_explicit_applicator::setup");

    n_domains = Fs.size();
    n_queues = (queues == nullptr ? 0 : queues->size());

    if(called_setup) eslog::error("setup was already called\n");
    if(apply_target != 'C' && apply_target != 'G') eslog::error("wrong apply target\n");
    if(vector_mem != 'C' && vector_mem != 'G') eslog::error("wrong vector memory\n");
    if(Fs_mem != 'C' && Fs_mem != 'G') eslog::error("wrong Fs memory\n");
    if((apply_target == 'G' || vector_mem == 'G') && !called_set_handles) eslog::error("handles must be set if gpu is used\n");
    if(D2C_old == nullptr) eslog::error("D2C is not set\n");
    if(D2C_old->size() != n_domains) eslog::error("wrong D2C\n");
    if(handles_dnblas != nullptr && handles_dnblas->size() != n_queues) eslog::error("wrong handles\n");
    if(n_domains == 0) eslog::error("no subdomains\n");
    if(omp_get_max_threads() > 1 && n_queues > 0 && omp_get_max_threads() != (int)n_queues) eslog::error("mismatch between num threads and num queues\n");
    if(Fs_mem == 'C' && std::any_of(Fs.begin(), Fs.end(), [](auto F){return !F->ator->is_data_accessible_cpu();}));
    if(Fs_mem == 'G' && std::any_of(Fs.begin(), Fs.end(), [](auto F){return !F->ator->is_data_accessible_gpu();}));
    if(!use_gpu && (Fs_mem == 'G' || vector_mem == 'G')) eslog::error("wrong memory, gpu is not used\n");
    if(!use_gpu && apply_target == 'G') eslog::error("cannot perform dual_operator application on GPU. GPU support not built or no GPU available.\n");

    wss_gpu_persistent = 0;

    ator_ws_gpu_persistent = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());
    ator_ws_gpu_tmp = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());

    need_copy_vectors = (apply_target != vector_mem);
    need_copy_Fs = (apply_target != Fs_mem);

    n_dofs_interfaces.resize(n_domains);
    total_dofs_interface = 0;
    for(size_t di = 0; di < n_domains; di++) {
        n_dofs_interfaces[di] += Fs[di]->nrows;
        total_dofs_interface += n_dofs_interfaces[di];
    }

    Allocator_new * ator_host = nullptr;
    if(apply_target == 'C') ator_host = AllocatorCPU_new::get_singleton();
    if(apply_target == 'G') ator_host = AllocatorHostPinned_new::get_singleton();

    Allocator_new * ator_target = nullptr;
    if(apply_target == 'C') ator_target = (use_gpu ? (Allocator_new*)AllocatorHostPinned_new::get_singleton() : (Allocator_new*)AllocatorCPU_new::get_singleton());
    if(apply_target == 'G') ator_target = ator_ws_gpu_persistent.get();

    D2C = MultiVectorDenseData_new<I,I>::convert_from(*D2C_old, ator_host);
    offsets = D2C.offsets;

    if(need_copy_Fs) {
        // todo: if hermitian and sharing allocation, copy and allocate only what i need
        Fs_2.resize(n_domains);
        for(size_t di = 0; di < n_domains; di++) {
            auto & F = *Fs[di];
            Fs_2[di].set(F.nrows, F.ncols, F.order, ator_target);
            Fs_2[di].prop = F.prop;
        }
    }

    xs_vec.set(n_domains, total_dofs_interface, ator_target);
    ys_vec.set(n_domains, total_dofs_interface, ator_target);

    if(apply_target == 'G') {
        d_D2C.set(D2C.num_vectors, D2C.size, ator_ws_gpu_persistent.get());
    }

    if(apply_target == 'G') {
        if(need_copy_Fs) {
            for(size_t di = 0; di < n_domains; di++) {
                wss_gpu_persistent += Fs_2[di].get_memory_impact();
            }
        }
        wss_gpu_persistent += xs_vec.get_memory_impact();
        wss_gpu_persistent += ys_vec.get_memory_impact();
        wss_gpu_persistent += d_D2C.get_memory_impact();
    }

    stacktimer::pop();
}



template<typename T, typename I>
size_t dualop_explicit_applicator<T,I>::get_wss_gpu_persistent()
{
    return wss_gpu_persistent;
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::set_ws_gpu_persistent(void * ws_gpu_persistent_)
{
    if(ws_gpu_persistent_ == nullptr && wss_gpu_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_gpu_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_gpu_persistent = ws_gpu_persistent_;
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::preprocess()
{
    stacktimer::push("dualop_explicit_applicator::preprocess");

    if(apply_target == 'G') {
        ator_ws_gpu_persistent->set(ws_gpu_persistent, wss_gpu_persistent);
    }

    if(need_copy_Fs) {
        for(size_t di = 0; di < n_domains; di++) {
            Fs_2[di].alloc();
        }
    }

    xs_vec.alloc();
    ys_vec.alloc();

    if(apply_target == 'G') {
        d_D2C.alloc();
        gpu::mgm::copy_submit(*main_q, D2C, d_D2C);
    }

    if(apply_target == 'C') {
        std::copy_n(offsets, n_domains + 1, xs_vec.offsets);
        std::copy_n(offsets, n_domains + 1, ys_vec.offsets);
    }
    if(apply_target == 'G') {
        gpu::mgm::copy_submit(*main_q, xs_vec.offsets, offsets, n_domains + 1);
        gpu::mgm::copy_submit(*main_q, ys_vec.offsets, xs_vec.offsets, n_domains + 1);
    }

    if(apply_target == 'G') {
        gpu::mgm::queue_wait(*main_q);
    }

    stacktimer::pop();
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::update_F(size_t di)
{
    if(need_copy_Fs) {
        gpu::mgm::copy_submit((*queues)[di % n_queues], *Fs[di], Fs_2[di]);
    }
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::apply(VectorDenseView_new<T> & x_cluster, VectorDenseView_new<T> & y_cluster, void * ws_gpu_tmp, size_t wss_gpu_tmp, const std::function<void(void)> & func_while_waiting)
{
    stacktimer::push("dualop_explicit_applicator::apply (vector)");

    if(x_cluster.size != y_cluster.size) eslog::error("size of x_cluster and y_cluster does not match\n");
    if(x_cluster.size != n_dofs_cluster_interface) eslog::error("incompatible cluster matrix size\n");

    if((x_cluster.ator->is_data_accessible_cpu() && vector_mem != 'C') || (x_cluster.ator->is_data_accessible_gpu() && vector_mem != 'G')) eslog::error("wrong x_cluster allocator\n");
    if((y_cluster.ator->is_data_accessible_cpu() && vector_mem != 'C') || (y_cluster.ator->is_data_accessible_gpu() && vector_mem != 'G')) eslog::error("wrong y_cluster allocator\n");

    stacktimer::push("dualop_apply_submit");

    ator_ws_gpu_tmp->set(ws_gpu_tmp, wss_gpu_tmp);

    Allocator_new * ator_tmp_target = nullptr;
    if(apply_target == 'C') ator_tmp_target = AllocatorHostPinned_new::get_singleton();
    if(apply_target == 'G') ator_tmp_target = ator_ws_gpu_tmp.get();

    VectorDenseData_new<T> x_cluster_2;
    VectorDenseData_new<T> y_cluster_2;
    if(need_copy_vectors) {
        x_cluster_2.set(x_cluster.size, ator_tmp_target);
        x_cluster_2.alloc();
        y_cluster_2.set(y_cluster.size, ator_tmp_target);
        y_cluster_2.alloc();
    }

    VectorDenseView_new<T> * x_cluster_to_use = (need_copy_vectors ? &x_cluster_2 : &x_cluster);
    VectorDenseView_new<T> * y_cluster_to_use = (need_copy_vectors ? &y_cluster_2 : &y_cluster);
    
    if(need_copy_vectors) {
        gpu::mgm::copy_submit(*main_q, x_cluster, x_cluster_2);
    }

    if(!timers_inner) stacktimer::disable();

    if(apply_target == 'C') {
        std::fill_n(y_cluster_to_use->vals, y_cluster_to_use->size, T{0});

        if(need_copy_Fs || need_copy_vectors) {
            gpu::mgm::device_wait(); // it would be better to wait for an event for each individual F, but whatever
        }

        #pragma omp parallel for schedule(static,1)
        for(size_t di = 0; di < n_domains; di++) {
            MatrixDenseView_new<T> * F_to_use = ((Fs_mem == 'C') ? Fs[di] : &Fs_2[di]);
            size_t n_dofs_interface = n_dofs_interfaces[di];

            for(size_t i = 0; i < n_dofs_interface; i++) {
                xs_vec.at(di,i) = x_cluster_to_use->vals[D2C.at(di,i)];
            }

            {
                VectorDenseView_new<T> x;
                VectorDenseView_new<T> y;
                x.set_view(n_dofs_interface, &xs_vec.at(di,0), xs_vec.ator);
                y.set_view(n_dofs_interface, &ys_vec.at(di,0), ys_vec.ator);

                if(is_hermitian<T>(F_to_use->prop.symm)) {
                    math::blas::apply_hermitian<T,I>(y, T{1}, *F_to_use, T{0}, x);
                }
                else {
                    math::blas::apply<T,I>(y, T{1}, *F_to_use, T{0}, x);
                }
            }

            for(size_t i = 0; i < n_dofs_interface; i++) {
                T val = ys_vec.at(di,i);
                T & dst = y_cluster_to_use->vals[D2C.at(di,i)];
                if constexpr(utils::is_complex<T>()) {
                    #pragma omp atomic
                    utils::real_ref(dst) += utils::real_ref(val);
                    #pragma omp atomic
                    utils::imag_ref(dst) += utils::imag_ref(val);
                }
                else {
                    #pragma omp atomic
                    dst += val;
                }
            }
        }
    }
    if(apply_target == 'G') {
        // scatter
        gpu::kernels::DCmap_scatter_new(*main_q, *x_cluster_to_use, xs_vec, d_D2C);

        gpu::mgm::queue_async_barrier({*main_q}, *queues);

        // apply
        #pragma omp parallel for schedule(static,1)
        for(size_t di = 0; di < n_domains; di++) {
            gpu::dnblas::handle & hd = (*handles_dnblas)[di % n_queues];

            MatrixDenseView_new<T> * F_to_use = ((Fs_mem == 'G') ? Fs[di] : &Fs_2[di]);
            size_t n_dofs_interface = n_dofs_interfaces[di];

            {
                T * x = xs_vec.vals + offsets[di];
                T * y = ys_vec.vals + offsets[di];

                if(is_hermitian<T>(F_to_use->prop.symm)) {
                    gpu::dnblas::hemv<T,I>(hd, n_dofs_interface, F_to_use->vals, F_to_use->ld, F_to_use->order, 'N', F_to_use->prop.uplo, x, y);
                }
                else {
                    gpu::dnblas::gemv<T,I>(hd, F_to_use->nrows, F_to_use->ncols, F_to_use->vals, F_to_use->ld, F_to_use->order, 'N', x, y);
                }
            }
        }

        // zerofill y_cluster on device
        gpu::mgm::memset_submit(*main_q, y_cluster_to_use->vals, y_cluster_to_use->size * sizeof(T), 0);

        gpu::mgm::queue_async_barrier(*queues, {*main_q});

        // gather
        gpu::kernels::DCmap_gather_new(*main_q, *y_cluster_to_use, ys_vec, d_D2C);
    }

    if(!timers_inner) stacktimer::enable();

    if(need_copy_vectors) {
        gpu::mgm::copy_submit(*main_q, y_cluster_2, y_cluster);
    }

    stacktimer::pop();

    stacktimer::push("dualop_apply_wait_intermediate");
    if(wait_intermediate && (apply_target == 'G' || need_copy_vectors)) {
        gpu::mgm::device_wait();
    }
    stacktimer::pop();

    stacktimer::push("dualop_apply_external_func");
    func_while_waiting();
    stacktimer::pop();

    stacktimer::push("dualop_apply_wait_final");
    if(apply_target == 'G' || need_copy_vectors) {
        gpu::mgm::device_wait();
    }
    stacktimer::pop();

    ator_ws_gpu_tmp->unset();

    stacktimer::pop();
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::apply(MatrixDenseView_new<T> & X_cluster, MatrixDenseView_new<T> & Y_cluster, void * ws_gpu_tmp, size_t wss_gpu_tmp, const std::function<void(void)> & func_while_waiting)
{
    stacktimer::push("dualop_explicit_applicator::apply (matrix)");

    if(X_cluster.nrows != Y_cluster.nrows || X_cluster.ncols != Y_cluster.ncols) eslog::error("size of X_cluster and Y_cluster does not match\n");
    if(X_cluster.order != Y_cluster.order) eslog::error("orders do not match\n");
    if(X_cluster.nrows != n_dofs_cluster_interface) eslog::error("incompatible cluster matrix size\n");

    if((X_cluster.ator->is_data_accessible_cpu() && vector_mem != 'C') || (X_cluster.ator->is_data_accessible_gpu() && vector_mem != 'G')) eslog::error("wrong X_cluster allocator\n");
    if((Y_cluster.ator->is_data_accessible_cpu() && vector_mem != 'C') || (Y_cluster.ator->is_data_accessible_gpu() && vector_mem != 'G')) eslog::error("wrong Y_cluster allocator\n");

    stacktimer::push("dualop_apply_submit");

    ator_ws_gpu_tmp->set(ws_gpu_tmp, wss_gpu_tmp);

    Allocator_new * ator_tmp_target = nullptr;
    if(apply_target == 'C') ator_tmp_target = AllocatorHostPinned_new::get_singleton();
    if(apply_target == 'G') ator_tmp_target = ator_ws_gpu_tmp.get();

    MatrixDenseData_new<T> X_cluster_2;
    MatrixDenseData_new<T> Y_cluster_2;
    if(need_copy_vectors) {
        X_cluster_2.set(X_cluster.nrows, X_cluster.ncols, X_cluster.order, ator_tmp_target);
        X_cluster_2.alloc();
        Y_cluster_2.set(Y_cluster.nrows, Y_cluster.ncols, Y_cluster.order, ator_tmp_target);
        Y_cluster_2.alloc();
    }

    MatrixDenseView_new<T> * X_cluster_to_use = (need_copy_vectors ? &X_cluster_2 : &X_cluster);
    MatrixDenseView_new<T> * Y_cluster_to_use = (need_copy_vectors ? &Y_cluster_2 : &Y_cluster);
    
    if(need_copy_vectors) {
        gpu::mgm::copy_submit(*main_q, X_cluster, X_cluster_2);
    }

    size_t cbmb_size = ator_ws_gpu_tmp->get_remaining_capacity();
    void * cbmb_mem = ator_ws_gpu_tmp->alloc(cbmb_size);
    AllocatorCBMB_new ator_gpu_cbmb(ator_ws_gpu_tmp.get(), cbmb_mem, cbmb_size);

    std::vector<MatrixDenseData_new<T>> Xs_gpu;
    std::vector<MatrixDenseData_new<T>> Ys_gpu;

    if(!timers_inner) stacktimer::disable();

    if(apply_target == 'C') {
        math::operations::fill_dnx<T>::do_all(Y_cluster_to_use, T{0});

        if(need_copy_Fs || need_copy_vectors) {
            gpu::mgm::device_wait(); // it would be better to wait for an event for each individual F, but whatever
        }

        #pragma omp parallel for schedule(static,1)
        for(size_t di = 0; di < n_domains; di++) {
            MatrixDenseView_new<T> * F_to_use = ((Fs_mem == 'C') ? Fs[di] : &Fs_2[di]);
            size_t n_dofs_interface = n_dofs_interfaces[di];

            MatrixDenseData_new<T> X;
            X.set(n_dofs_interface, X_cluster_to_use->ncols, X_cluster_to_use->order, AllocatorCPU_new::get_singleton());
            X.alloc();
            
            MatrixDenseData_new<T> Y;
            Y.set(n_dofs_interface, Y_cluster_to_use->ncols, Y_cluster_to_use->order, AllocatorCPU_new::get_singleton());
            Y.alloc();

            VectorDenseView_new<I> my_D2C;
            my_D2C.set_view(n_dofs_interface, &D2C.at(di,0), D2C.ator);

            math::operations::submatrix_dnx_dnx_noncontig<T,I>::do_all(X_cluster_to_use, &X, &my_D2C, nullptr);

            if(is_hermitian<T>(F_to_use->prop.symm)) {
                math::blas::hemm<T>(*F_to_use, X, Y);
            }
            else {
                math::blas::gemm<T>(*F_to_use, X, Y);
            }

            math::operations::supermatrix_dnx_dnx_noncontig<T,I>::do_all(&Y, Y_cluster_to_use, &my_D2C, nullptr, math::operations::supermatrix_dnx_dnx_noncontig<T,I>::mode::accumulate_atomic);
        }
    }
    if(apply_target == 'G') {
        gpu::mgm::memset_submit(*main_q, Y_cluster_to_use->vals, Y_cluster_to_use->get_size_primary() * Y_cluster_to_use->ld * sizeof(T), 0);
        gpu::mgm::queue_async_barrier({*main_q}, *queues);

        Xs_gpu.resize(n_domains);
        Ys_gpu.resize(n_domains);

        #pragma omp parallel for schedule(static,1)
        for(size_t di = 0; di < n_domains; di++) {
            gpu::mgm::queue & q = (*queues)[di % n_queues];
            gpu::dnblas::handle & hd = (*handles_dnblas)[di % n_queues];

            MatrixDenseView_new<T> * F_to_use = ((Fs_mem == 'G') ? Fs[di] : &Fs_2[di]);
            size_t n_dofs_interface = n_dofs_interfaces[di];

            MatrixDenseData_new<T> & X = Xs_gpu[di];
            MatrixDenseData_new<T> & Y = Ys_gpu[di];

            X.set(n_dofs_interface, X_cluster_to_use->ncols, X_cluster_to_use->order, &ator_gpu_cbmb);
            Y.set(n_dofs_interface, Y_cluster_to_use->ncols, Y_cluster_to_use->order, &ator_gpu_cbmb);

            std::unique_ptr<gpu::operations::hemm_ddnx_ddny_ddnz<T>> op_hemm;
            std::unique_ptr<gpu::operations::gemm_ddnx_ddny_ddnz<T>> op_gemm;

            size_t wss_tmp_ghemm;
            if(is_hermitian<T>(F_to_use->prop.symm)) {
                op_hemm = gpu::operations::hemm_ddnx_ddny_ddnz<T>::make();
                op_hemm->set_handles(q, hd);
                op_hemm->set_matrix_A(F_to_use);
                op_hemm->set_matrix_B(&X);
                op_hemm->set_matrix_C(&Y);
                op_hemm->set_coefficients(T{1}, T{0});
                op_hemm->setup();
                wss_tmp_ghemm = op_hemm->get_wss_tmp_perform();
            }
            else {
                op_gemm = gpu::operations::gemm_ddnx_ddny_ddnz<T>::make();
                op_gemm->set_handles(q, hd);
                op_gemm->set_matrix_A(F_to_use);
                op_gemm->set_matrix_B(&X);
                op_gemm->set_matrix_C(&Y);
                op_gemm->set_coefficients(T{1}, T{0});
                op_gemm->setup();
                wss_tmp_ghemm = op_gemm->get_wss_tmp_perform();
            }

            void * ws_tmp_ghemm = nullptr;
            ator_gpu_cbmb.resource.do_transaction([&](){
                X.alloc();
                Y.alloc();
                ws_tmp_ghemm = ator_gpu_cbmb.alloc(wss_tmp_ghemm);
            });

            VectorDenseView_new<I> d_my_D2C;
            d_my_D2C.set_view(n_dofs_interface, d_D2C.vals + offsets[di], d_D2C.ator);

            gpu::operations::submatrix_ddnx_ddnx_noncontig<T,I>::submit_all(q, X_cluster_to_use, &X, &d_my_D2C, nullptr);

            if(is_hermitian<T>(F_to_use->prop.symm)) {
                op_hemm->perform_submit(ws_tmp_ghemm);
            }
            else {
                op_gemm->perform_submit(ws_tmp_ghemm);
            }

            gpu::operations::supermatrix_ddnx_ddnx_noncontig<T,I>::submit_all(q, &Y, Y_cluster_to_use, &d_my_D2C, nullptr, gpu::operations::supermatrix_ddnx_ddnx_noncontig<T,I>::mode::accumulate_atomic);

            gpu::mgm::submit_host_function(q, [&,ws_tmp_ghemm,di](){
                Xs_gpu[di].free();
                Ys_gpu[di].free();
                void * ws_tmp_ghemm_ = ws_tmp_ghemm;
                ator_gpu_cbmb.free(ws_tmp_ghemm_);
            });
        }
    }

    if(!timers_inner) stacktimer::enable();

    gpu::mgm::queue_async_barrier(*queues, {*main_q});

    if(need_copy_vectors) {
        gpu::mgm::copy_submit(*main_q, Y_cluster_2, Y_cluster);
    }

    stacktimer::pop();

    stacktimer::push("dualop_apply_wait_intermediate");
    if(wait_intermediate && (apply_target == 'G' || need_copy_vectors)) {
        gpu::mgm::device_wait();
    }
    stacktimer::pop();

    stacktimer::push("dualop_apply_external_func");
    func_while_waiting();
    stacktimer::pop();

    stacktimer::push("dualop_apply_wait_final");
    if(apply_target == 'G' || need_copy_vectors) {
        gpu::mgm::device_wait();
    }
    stacktimer::pop();

    ator_ws_gpu_tmp->unset();

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class dualop_explicit_applicator<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        INSTANTIATE_T(std::complex<double>)

    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
