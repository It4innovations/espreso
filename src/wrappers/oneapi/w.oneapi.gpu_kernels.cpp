
#ifdef HAVE_ONEAPI
#ifdef ESPRESO_USE_WRAPPER_GPU_ONEAPI

#include "gpu/gpu_kernels.h"
#include "w.oneapi.gpu_management.h"
#include "basis/utilities/utils.h"

namespace espreso {
namespace gpu {
namespace kernels {

    namespace
    {
        template<typename T>
        void my_atomicadd(T * dst, T val)
        {
            sycl::atomic_ref<T, sycl::memory_order::relaxed, sycl::memory_scope::device, sycl::access::address_space::global_space> dst_ref(*dst);
            dst_ref += val;
        }
        template<typename T>
        void my_atomicadd(std::complex<T> * dst, std::complex<T> val)
        {
            my_atomicadd(&utils::real_ref(*dst), utils::real_ref(val));
            my_atomicadd(&utils::imag_ref(*dst), utils::imag_ref(val));
        }
    }

    template<typename T, typename I, typename A>
    void DCmap_scatter(mgm::queue & q, Vector_Dense<T*,I,A> & domain_vector_pointers, const Vector_Dense<I,I,A> & n_dofs_interfaces, const Vector_Dense<T,I,A> & cluster_vector, const Vector_Dense<I*,I,A> & D2Cs)
    {
        static_assert(A::is_data_device_accessible, "data has to be device accessible");
        int wpw = 256;
        I n_domains = domain_vector_pointers.size;
        sycl::nd_range<1> range(sycl::range<1>(n_domains * wpw), sycl::range<1>(wpw));
        T ** domain_vectors = domain_vector_pointers.vals;
        I * n_dofs_interfaces_vals = n_dofs_interfaces.vals;
        T * cluster_vector_vals = cluster_vector.vals;
        I ** D2Cs_vals = D2Cs.vals;
        q->q.parallel_for(
            range,
            [=](sycl::nd_item<1> item) {
                sycl::group g = item.get_group();
                I di = g.get_group_linear_id();
                I n_dofs_interface = n_dofs_interfaces_vals[di];
                T * domain_vector = domain_vectors[di];
                I * D2C = D2Cs_vals[di];
                for(I i = g.get_local_linear_id(); i < n_dofs_interface; i += g.get_local_linear_range()) {
                    domain_vector[i] = cluster_vector_vals[D2C[i]];
                }
            }
        );
    }

    template<typename T, typename I>
    void DCmap_scatter_new(mgm::queue & q, const VectorDenseView_new<T> & vec_cluster, MultiVectorDenseView_new<T,I> & vecs_subdomains, const MultiVectorDenseView_new<I,I> & D2C)
    {
        if(!vec_cluster.ator->is_data_accessible_gpu()) eslog::error("wrong allocator\n");
        if(!vecs_subdomains.ator->is_data_accessible_gpu()) eslog::error("wrong allocator\n");
        if(!D2C.ator->is_data_accessible_gpu()) eslog::error("wrong allocator\n");
        int wpw = 256;
        I n_domains = vecs_subdomains.num_vectors;
        sycl::nd_range<1> range(sycl::range<1>(n_domains * wpw), sycl::range<1>(wpw));
        const T * cluster_vector_vals = vec_cluster.vals;
        T * domain_vectors_vals = vecs_subdomains.vals;
        const I * domain_vectors_offsets = vecs_subdomains.offsets;
        const I * D2C_vals = D2C.vals;
        const I * D2C_offsets = D2C.offsets;
        q->q.parallel_for(
            range,
            [=](sycl::nd_item<1> item) {
                sycl::group g = item.get_group();
                I di = g.get_group_linear_id();
                I start = domain_vectors_offsets[di];
                I end = domain_vectors_offsets[di+1];
                I size = end - start;
                T * vec = domain_vectors_vals + start;
                const I * D2C_domain = D2C_vals + D2C_offsets[di];
                for(I i = g.get_local_linear_id(); i < size; i += g.get_local_linear_range()) {
                    vec[i] = cluster_vector_vals[D2C_domain[i]];
                }
            }
        );
    }

    template<typename T, typename I, typename A>
    void DCmap_gather(mgm::queue & q, const Vector_Dense<T*,I,A> & domain_vector_pointers, const Vector_Dense<I,I,A> & n_dofs_interfaces, Vector_Dense<T,I,A> & cluster_vector, const Vector_Dense<I*,I,A> & D2Cs)
    {
        static_assert(A::is_data_device_accessible, "data has to be device accessible");
        int wpw = 256;
        I n_domains = domain_vector_pointers.size;
        sycl::nd_range<1> range(sycl::range<1>(n_domains * wpw), sycl::range<1>(wpw));
        T ** domain_vectors = domain_vector_pointers.vals;
        I * n_dofs_interfaces_vals = n_dofs_interfaces.vals;
        T * cluster_vector_vals = cluster_vector.vals;
        I ** D2Cs_vals = D2Cs.vals;
        q->q.parallel_for(
            range,
            [=](sycl::nd_item<1> item) {
                sycl::group g = item.get_group();
                I di = g.get_group_linear_id();
                I n_dofs_interface = n_dofs_interfaces_vals[di];
                T * domain_vector = domain_vectors[di];
                I * D2C = D2Cs_vals[di];
                for(I i = g.get_local_linear_id(); i < n_dofs_interface; i += g.get_local_linear_range()) {
                    my_atomicadd(&cluster_vector_vals[D2C[i]], domain_vector[i]);
                }
            }
        );
    }

    template<typename T, typename I>
    void DCmap_gather_new(mgm::queue & q, VectorDenseView_new<T> & vec_cluster, const MultiVectorDenseView_new<T,I> & vecs_subdomains, const MultiVectorDenseView_new<I,I> & D2C)
    {
        if(!vec_cluster.ator->is_data_accessible_gpu()) eslog::error("wrong allocator\n");
        if(!vecs_subdomains.ator->is_data_accessible_gpu()) eslog::error("wrong allocator\n");
        if(!D2C.ator->is_data_accessible_gpu()) eslog::error("wrong allocator\n");
        int wpw = 256;
        I n_domains = vecs_subdomains.num_vectors;
        sycl::nd_range<1> range(sycl::range<1>(n_domains * wpw), sycl::range<1>(wpw));
        T * cluster_vector_vals = vec_cluster.vals;
        const T * domain_vectors_vals = vecs_subdomains.vals;
        const I * domain_vectors_offsets = vecs_subdomains.offsets;
        const I * D2C_vals = D2C.vals;
        const I * D2C_offsets = D2C.offsets;
        q->q.parallel_for(
            range,
            [=](sycl::nd_item<1> item) {
                sycl::group g = item.get_group();
                I di = g.get_group_linear_id();
                I start = domain_vectors_offsets[di];
                I end = domain_vectors_offsets[di+1];
                I size = end - start;
                const T * vec = domain_vectors_vals + start;
                const I * D2C_domain = D2C_vals + D2C_offsets[di];
                for(I i = g.get_local_linear_id(); i < size; i += g.get_local_linear_range()) {
                    my_atomicadd(&cluster_vector_vals[D2C_domain[i]], vec[i]);
                }
            }
        );
    }

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_matrix_triangle(mgm::queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input, char fill, char order)
    {
        static_assert(Ao::is_data_device_accessible, "matrix data has to be device accessible");
        static_assert(Ai::is_data_device_accessible, "matrix data has to be device accessible");
        if(output.nrows != input.nrows || output.ncols != input.ncols || input.nrows != input.ncols) eslog::error("matrix dimensions do not match\n");

        if(order == 'R') {
            T * output_vals = output.vals;
            T * input_vals = input.vals;
            I output_ld = output.get_ld();
            I input_ld = input.get_ld();
            I n = input.nrows;
            int wpw = 256;
            int ngroups = (n - 1) / 2 + 1;
            q->q.parallel_for(
                sycl::nd_range(sycl::range<1>(ngroups * wpw), sycl::range<1>(wpw)),
                [=](sycl::nd_item<1> item) {
                    // assuming row-major matrix
                    // launch with (n-1)/2+1 blocks (one block per two rows)
                    // each block handles rows blockIdx.x (for stage=0) and n-blockIdx.x-1 (for stage=1).
                    // in stage=0, warps are normally sequential next to each other
                    // in stage=1, the order of warps is reversed for better load balance between warps

                    sycl::sub_group sg = item.get_sub_group();
                    I wgidx = item.get_group_linear_id();
                    I wgsize = item.get_local_range().size();
                    I sgidx = sg.get_group_linear_id();
                    I sgcount = sg.get_group_linear_range();
                    I sgsize = sg.get_max_local_range()[0];
                    I idxinsg = sg.get_local_linear_id();
                    #pragma unroll
                    for(int stage = 0; stage <= 1; stage++) {
                        if(stage == 1 && wgidx == n/2) return;
                        I r = (stage == 0) ? (wgidx) : (n - wgidx - 1);
                        I chunk_idx = (stage == 0) ? (sgidx) : (sgcount - sgidx - 1);
                        I myidx = chunk_idx * sgsize + idxinsg;
                        const T * row_in = input_vals + r * input_ld;
                        T * row_out = output_vals + r * output_ld;
                        if(fill == 'L') {
                            for(I c = myidx; c <= r; c += wgsize) row_out[c] = row_in[c];
                        }
                        if(fill == 'U') {
                            I skip_leaps = r / wgsize;
                            for(I c = skip_leaps * wgsize + myidx; c < n; c += wgsize) if(c >= r) row_out[c] = row_in[c];
                        }
                    }
                    
                }
            );
        }
        else if(order == 'C') {
            char fill_compl = mgm::fill_change(fill);
            copy_matrix_triangle(q, output, input, fill_compl, 'R');
        }
        else {
            eslog::error("invalid order %c\n", order);
        }
    }

}
}
}

#include "gpu/gpu_kernels.inst.hpp"

#endif
#endif
