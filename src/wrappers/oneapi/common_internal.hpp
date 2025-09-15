
#ifndef SRC_WRAPPERS_ONEAPI_COMMON_INTERNAL_HPP
#define SRC_WRAPPERS_ONEAPI_COMMON_INTERNAL_HPP



namespace espreso {
namespace gpu {

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
}



#endif /* SRC_WRAPPERS_ONEAPI_COMMON_INTERNAL_HPP */
