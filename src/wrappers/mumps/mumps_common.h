
#ifndef SRC_WRAPPERS_MUMPS_MUMPS_COMMON_H
#define SRC_WRAPPERS_MUMPS_MUMPS_COMMON_H

#include <smumps_c.h>
#include <dmumps_c.h>
#include <cmumps_c.h>
#include <zmumps_c.h>

#include <complex>
#include "math/primitives_new.h"



namespace espreso {



template<typename T> struct my_mumps_handle;

template<> struct my_mumps_handle<float> { using type = SMUMPS_STRUC_C; };
template<> struct my_mumps_handle<double> { using type = DMUMPS_STRUC_C; };
template<> struct my_mumps_handle<std::complex<float>> { using type = CMUMPS_STRUC_C; };
template<> struct my_mumps_handle<std::complex<double>> { using type = ZMUMPS_STRUC_C; };

template<typename T> using my_mumps_handle_t = typename my_mumps_handle<T>::type;



template<typename T> struct cpp_type_to_mumps_type;

template<> struct cpp_type_to_mumps_type<float> { using type = float; };
template<> struct cpp_type_to_mumps_type<double> { using type = double; };
template<> struct cpp_type_to_mumps_type<std::complex<float>> { using type = mumps_complex; };
template<> struct cpp_type_to_mumps_type<std::complex<double>> { using type = mumps_double_complex; };

template<typename T> using cpp_type_to_mumps_type_t = typename cpp_type_to_mumps_type<T>::type;



template<typename T>
static inline void call_mumps(my_mumps_handle_t<T> & handle)
{
    // mumps is not thread safe. huge performance hit, but at least it works
    #pragma omp critical(espreso_call_mumps)
    {
        if constexpr(std::is_same_v<T,float>) smumps_c(&handle);
        if constexpr(std::is_same_v<T,double>) dmumps_c(&handle);
        if constexpr(std::is_same_v<T,std::complex<float>>) cmumps_c(&handle);
        if constexpr(std::is_same_v<T,std::complex<double>>) zmumps_c(&handle);
        
        int status = handle.info[0];
        if(status != 0) {
            int detail = handle.info[1];
            eslog::error("mumps error %d, detail %d\n", status, detail);
        }
    }
}



template<typename T, typename I>
void mumps_helper_csx_to_ijv(MatrixCsxView_new<T,I> & M, VectorDenseData_new<I> & ijv_rowidxs, VectorDenseData_new<I> & ijv_colidxs);



}



#endif /* SRC_WRAPPERS_MUMPS_MUMPS_COMMON_H */
