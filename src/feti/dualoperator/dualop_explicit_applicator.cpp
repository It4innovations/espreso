
#include "dualop_explicit_applicator.h"
#include "math/primitives_new/allocator_new.h"
#include "gpu/gpu_kernels.h"
#include "basis/utilities/stacktimer.h"

#include <numeric>



namespace espreso {



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::set_config(bool parallel_apply_)
{
    parallel_apply = parallel_apply_;
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
void dualop_explicit_applicator<T,I>::set_dimensions(FETI<T> & feti)
{
    n_dofs_cluster_interface = feti.lambdas.size;
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::set_vector_memory(char vector_mem_)
{
    if(called_setup) eslog::error("cannot change vector memory anymore\n");

    vector_mem = vector_mem_;
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
    n_domains = Fs.size();
    n_queues = (queues == nullptr ? 0 : queues->size());

    if(called_setup) eslog::error("setup was already called\n");
    if(apply_target != 'C' && apply_target != 'G') eslog::error("wrong apply target\n");
    if(vector_mem != 'C' && vector_mem != 'G') eslog::error("wrong vector memory\n");
    if((apply_target == 'G' || vector_mem == 'G') && !called_set_handles) eslog::error("handles must be set if gpu is used\n");
    if(D2C_old == nullptr) eslog::error("D2C is not set\n");
    if(D2C_old->size() != n_domains) eslog::error("wrong D2C\n");
    if(handles_dnblas != nullptr && handles_dnblas->size() != n_queues) eslog::error("wrong handles\n");
    if(n_domains == 0) eslog::error("no subdomains\n");
    if(omp_get_max_threads() > 1 && n_queues > 0 && omp_get_max_threads() != (int)n_queues) eslog::error("mismatch between num threads and num queues\n");

    if(Fs[0]->ator->is_data_accessible_cpu()) Fs_mem = 'C';
    if(Fs[0]->ator->is_data_accessible_gpu()) Fs_mem = 'G';
    for(auto F : Fs) if(F->ator->is_data_accessible_cpu() != Fs[0]->ator->is_data_accessible_cpu() || F->ator->is_data_accessible_gpu() != Fs[0]->ator->is_data_accessible_gpu()) eslog::error("different allocators for different F\n");

    need_copy_vectors = (apply_target != vector_mem);
    need_copy_Fs = (apply_target != Fs_mem);

    n_dofs_interfaces.resize(n_domains);
    total_dofs_interface = 0;
    for(size_t di = 0; di < n_domains; di++) {
        n_dofs_interfaces[di] += Fs[di]->nrows;
        total_dofs_interface += n_dofs_interfaces[di];
    }

    Allocator_new * ator_2 = nullptr;
    if(apply_target == 'C') ator_2 = AllocatorCPU_new::get_singleton();
    if(apply_target == 'G') ator_2 = AllocatorGPU_new::get_singleton();

    Allocator_new * ator_host = nullptr;
    if(apply_target == 'C') ator_host = AllocatorCPU_new::get_singleton();
    if(apply_target == 'G') ator_host = AllocatorHostPinned_new::get_singleton();

    D2C = MultiVectorDenseData_new<I,I>::convert_from(*D2C_old, ator_host);
    offsets = D2C.offsets;

    if(need_copy_Fs) {
        // todo: if hermitian and sharing allocation, copy and allocate only what i need
        Fs_2.resize(n_domains);
        for(size_t di = 0; di < n_domains; di++) {
            auto & F = *Fs[di];
            Fs_2[di].set(F.nrows, F.ncols, F.order, ator_2);
            Fs_2[di].prop = F.prop;
            Fs_2[di].conj = F.conj;
        }
    }
    
    if(need_copy_vectors) {
        x_cluster_2.set(n_dofs_cluster_interface, ator_2);
        y_cluster_2.set(n_dofs_cluster_interface, ator_2);
    }

    xs.set(n_domains, total_dofs_interface, ator_2);
    ys.set(n_domains, total_dofs_interface, ator_2);

    if(apply_target == 'G') {
        d_D2C.set(D2C.num_vectors, D2C.size, AllocatorGPU_new::get_singleton());
    }

    if(apply_target == 'G') {
        if(need_copy_Fs) {
            for(size_t di = 0; di < n_domains; di++) {
                gpu_wss_internal += Fs_2[di].get_memory_impact();
            }
        }
        if(need_copy_vectors) {
            gpu_wss_internal += x_cluster_2.get_memory_impact();
            gpu_wss_internal += y_cluster_2.get_memory_impact();
        }
        gpu_wss_internal += xs.get_memory_impact();
        gpu_wss_internal += ys.get_memory_impact();
        gpu_wss_internal += d_D2C.get_memory_impact();
    }
}



template<typename T, typename I>
size_t dualop_explicit_applicator<T,I>::get_gpu_wss_internal()
{
    return gpu_wss_internal;
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::preprocess()
{
    if(need_copy_Fs) {
        for(size_t di = 0; di < n_domains; di++) {
            Fs_2[di].alloc();
        }
    }

    if(need_copy_vectors) {
        x_cluster_2.alloc();
        y_cluster_2.alloc();
    }

    xs.alloc();
    ys.alloc();

    if(apply_target == 'G') {
        d_D2C.alloc();
        gpu::mgm::copy_submit(*main_q, D2C, d_D2C);
    }

    if(apply_target == 'C') {
        std::copy_n(offsets, n_domains + 1, xs.offsets);
        std::copy_n(offsets, n_domains + 1, ys.offsets);
    }
    if(apply_target == 'G') {
        gpu::mgm::copy_submit(*main_q, xs.offsets, offsets, n_domains + 1);
        gpu::mgm::copy_submit(*main_q, ys.offsets, xs.offsets, n_domains + 1);
    }

    if(apply_target == 'G') {
        gpu::mgm::queue_wait(*main_q);
    }
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::update_F(size_t di)
{
    if(need_copy_Fs) {
        gpu::mgm::copy_submit((*queues)[di % n_queues], *Fs[di], Fs_2[di]);
    }
}



template<typename T, typename I>
void dualop_explicit_applicator<T,I>::apply(VectorDenseView_new<T> & x_cluster, VectorDenseView_new<T> & y_cluster)
{
    stacktimer::push("dualop_explicit_applicator::apply");

    if((x_cluster.ator->is_data_accessible_cpu() && vector_mem != 'C') || (x_cluster.ator->is_data_accessible_gpu() && vector_mem != 'G')) eslog::error("wrong x_cluster allocator\n");
    if((y_cluster.ator->is_data_accessible_cpu() && vector_mem != 'C') || (y_cluster.ator->is_data_accessible_gpu() && vector_mem != 'G')) eslog::error("wrong y_cluster allocator\n");

    VectorDenseView_new<T> * x_cluster_to_use = (need_copy_vectors ? &x_cluster_2 : &x_cluster);
    VectorDenseView_new<T> * y_cluster_to_use = (need_copy_vectors ? &y_cluster_2 : &y_cluster);
    
    if(need_copy_vectors) {
        gpu::mgm::copy_submit(*main_q, x_cluster, x_cluster_2);
    }

    if(apply_target == 'C') {
        std::fill_n(y_cluster_to_use->vals, y_cluster_to_use->size, T{0});

        if(need_copy_Fs || need_copy_vectors) {
            gpu::mgm::device_wait(); // it would be better to wait for an event for each individual F, but whatever
        }

        #pragma omp parallel for schedule(static,1) if(parallel_apply)
        for(size_t di = 0; di < n_domains; di++) {
            MatrixDenseView_new<T> * F_to_use = ((Fs_mem == 'C') ? Fs[di] : &Fs_2[di]);
            size_t n_dofs_interface = n_dofs_interfaces[di];

            for(size_t i = 0; i < n_dofs_interface; i++) {
                xs.at(di,i) = x_cluster_to_use->vals[D2C.at(di,i)];
            }

            {
                VectorDenseView_new<T> x;
                VectorDenseView_new<T> y;
                x.set_view(n_dofs_interface, &xs.at(di,0), xs.ator);
                y.set_view(n_dofs_interface, &ys.at(di,0), ys.ator);

                if(is_hermitian<T>(F_to_use->prop.symm)) {
                    math::blas::apply_hermitian<T,I>(y, T{1}, *F_to_use, T{0}, x);
                }
                else {
                    math::blas::apply<T,I>(y, T{1}, *F_to_use, T{0}, x);
                }
            }

            for(size_t i = 0; i < n_dofs_interface; i++) {
                T val = ys.at(di,i);
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
        gpu::kernels::DCmap_scatter_new(*main_q, *x_cluster_to_use, xs, d_D2C);

        gpu::mgm::queue_async_barrier({*main_q}, *queues);

        // apply
        #pragma omp parallel for schedule(static,1) if(parallel_apply)
        for(size_t di = 0; di < n_domains; di++) {
            gpu::dnblas::handle & hd = (*handles_dnblas)[di % n_queues];

            MatrixDenseView_new<T> * F_to_use = ((Fs_mem == 'G') ? Fs[di] : &Fs_2[di]);
            size_t n_dofs_interface = n_dofs_interfaces[di];

            {
                T * x = xs.vals + offsets[di];
                T * y = ys.vals + offsets[di];

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
        gpu::kernels::DCmap_gather_new(*main_q, *y_cluster_to_use, ys, d_D2C);
    }

    if(need_copy_vectors) {
        gpu::mgm::copy_submit(*main_q, y_cluster_2, y_cluster);
    }

    if(apply_target == 'G' || need_copy_vectors) {
        gpu::mgm::device_wait();
    }

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
