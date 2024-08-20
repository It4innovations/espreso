
#ifdef HAVE_ONEAPI

#include "gpu/gpu_kernels.h"

namespace espreso {
namespace gpu {
namespace kernels {

    namespace
    {
        template<typename T>
        void my_atomicadd(T * dst, T val)
        {
            sycl::atomic_ref<T, sycl::memory_order::relaxed, sycl::memory_scope::device, sycl::access::address_space::global_space> dst_ref(&val);
            dst_ref += val;
        }
        template<typename T>
        void my_atomicadd(std::complex<T> * dst, std::complex<T> val)
        {
            my_atomicadd(&utils::real_ref(*dst), utils::real_ref(val));
            my_atomicadd(&utils::imag_ref(*dst), utils::imag_ref(val));
        }
    }

    template<typename T, typename I>
    void DCmap_scatter(mgm::queue & q, Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, const Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs)
    {
        int wpw = 256;
        I n_domains = domain_vector_pointers.size;
        sycl::nd_range<1> range(sycl::range<1>(n_domains * wpw), sycl::range<1>(wpw));
        T ** domain_vectors = domain_vector_pointers.vals;
        I * n_dofs_interfaces_vals = n_dofs_interfaces.vals;
        T * cluster_vector_vals = cluster_vector.vals;
        I * D2Cs_vals = D2Cs.vals;
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
    void DCmap_gather(mgm::queue & q, const Vector_Dense<T*,I,mgm::Ad> & domain_vector_pointers, const Vector_Dense<I,I,mgm::Ad> & n_dofs_interfaces, Vector_Dense<T,I,mgm::Ad> & cluster_vector, const Vector_Dense<I*,I,mgm::Ad> & D2Cs)
    {
        int wpw = 256;
        I n_domains = domain_vector_pointers.size;
        sycl::nd_range<1> range(sycl::range<1>(n_domains * wpw), sycl::range<1>(wpw));
        T ** domain_vectors = domain_vector_pointers.vals;
        I * n_dofs_interfaces_vals = n_dofs_interfaces.vals;
        T * cluster_vector_vals = cluster_vector.vals;
        I * D2Cs_vals = D2Cs.vals;
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

    template<typename T, typename I, typename Ao, typename Ai>
    void copy_matrix_triangle(mgm::queue & q, Matrix_Dense<T,I,Ao> & output, const Matrix_Dense<T,I,Ai> & input, char fill, char order)
    {

    }

}
}
}

#include "gpu/gpu_kernels.inst.hpp"

#endif
