
namespace espreso {
namespace gpu {
namespace mgm {

    #define INSTANTIATE_T_I_AHOST_ADEVICE(T,I,Ahost,Adevice) \
    template void copy_submit_h2d<T,I,Adevice,Ahost  >(queue & q, Vector_Dense<T,I,Adevice> & output, const Vector_Dense<T,I,Ahost>   & input); \
    template void copy_submit_d2h<T,I,Ahost,  Adevice>(queue & q, Vector_Dense<T,I,Ahost>   & output, const Vector_Dense<T,I,Adevice> & input); \
    template void copy_submit_h2d<T,I,Adevice,Ahost  >(queue & q, Matrix_Dense<T,I,Adevice> & output, const Matrix_Dense<T,I,Ahost>   & input); \
    template void copy_submit_d2h<T,I,Ahost,  Adevice>(queue & q, Matrix_Dense<T,I,Ahost>   & output, const Matrix_Dense<T,I,Adevice> & input); \
    template void copy_submit_h2d<T,I,Adevice,Ahost  >(queue & q, Matrix_CSR<T,I,Adevice> & output, const Matrix_CSR<T,I,Ahost>   & input, bool copy_pattern, bool copy_vals); \
    template void copy_submit_d2h<T,I,Ahost,  Adevice>(queue & q, Matrix_CSR<T,I,Ahost>   & output, const Matrix_CSR<T,I,Adevice> & input, bool copy_pattern, bool copy_vals);

        #define INSTANTIATE_T_I(T,I) \
        INSTANTIATE_T_I_AHOST_ADEVICE(T, I, mgm::Ah,       mgm::Ad) \
        INSTANTIATE_T_I_AHOST_ADEVICE(T, I, mgm::Ah,       cbmba_d) \
        INSTANTIATE_T_I_AHOST_ADEVICE(T, I, cpu_allocator, mgm::Ad) \
        INSTANTIATE_T_I_AHOST_ADEVICE(T, I, cpu_allocator, cbmba_d)

            #define INSTANTIATE_T(T) \
            template void copy_submit_h2d<T>(queue & q, T * dst, T const * src, size_t num_elements); \
            template void copy_submit_d2h<T>(queue & q, T * dst, T const * src, size_t num_elements); \
            INSTANTIATE_T_I(T, int32_t) \
            /* INSTANTIATE_T_I(T, int64_t) */

                // INSTANTIATE_T(float)
                INSTANTIATE_T(double)
                // INSTANTIATE_T(std::complex<float>)
                // INSTANTIATE_T(std::complex<double>)
                // INSTANTIATE_T(float*)
                INSTANTIATE_T(double*)
                // INSTANTIATE_T(std::complex<float>*)
                // INSTANTIATE_T(std::complex<double>*)
                INSTANTIATE_T(int32_t)
                INSTANTIATE_T(int32_t*)
                // INSTANTIATE_T(int64_t)
                // INSTANTIATE_T(int64_t*)

            #undef INSTANTIATE_T
        #undef INSTANTIATE_T_I
    #undef INSTANTIATE_T_I_AHOST_ADEVICE

}
}
}
