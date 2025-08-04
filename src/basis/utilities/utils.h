
#ifndef BASIS_UTILITIES_UTILS_H_
#define BASIS_UTILITIES_UTILS_H_

#include <cstddef>
#include <vector>
#include <limits>
#include <complex>
#include <algorithm>

namespace espreso {
namespace utils {

    template<typename Ttype>
    Ttype sizesToOffsets(Ttype *begin, Ttype *end, Ttype offset = 0);

    template<typename Ttype>
    Ttype sizesToOffsets(std::vector<Ttype> &sizes, Ttype offset = 0)
    {
        return sizesToOffsets(sizes.data(), sizes.data() + sizes.size(), offset);
    }

    template<typename Ttype>
    std::vector<Ttype> sizesToOffsets(std::vector<std::vector<Ttype> > &sizes, const std::vector<Ttype> &offsets);

    template<typename Ttype>
    std::vector<Ttype> sizesToOffsets(std::vector<std::vector<Ttype> > &sizes)
    {
        return sizesToOffsets(sizes, std::vector<Ttype>(sizes.front().size()));
    }

    template<typename Ttype, typename Tpermutation>
    void permute(std::vector<Ttype> &data, const std::vector<Tpermutation> &permutation, size_t elementsize = 1);

    template<typename Ttype>
    void threadDistributionToFullDistribution(std::vector<std::vector<Ttype> > &distribution);

    template<typename Ttype>
    void threadDistributionToFullDistribution(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

    template<typename Ttype>
    void removeDuplicates(std::vector<Ttype> &data, size_t begin = 0);

    template<typename Ttype>
    void sortAndRemoveDuplicates(std::vector<Ttype> &data, size_t begin = 0);

    template<typename Ttype>
    void sortAndRemoveDuplicates(std::vector<std::vector<Ttype> > &data);

    template<typename Ttype>
    void mergeThreadedUniqueData(std::vector<std::vector<Ttype> > &data);

    template<typename Ttype>
    void mergeThreadedUniqueData(std::vector<std::vector<std::vector<Ttype> > > &data);

    template<typename Ttype>
    void inplaceMerge(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

    template<typename Ttype>
    void inplaceMerge(std::vector<std::vector<Ttype> > &data);

    template<typename Ttype>
    void sortWithInplaceMerge(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

    template<typename Ttype>
    void sortWithInplaceMerge(std::vector<std::vector<Ttype> > &data);

    template<typename Ttype>
    void sortWithUniqueMerge(std::vector<std::vector<Ttype> > &data);

    template<typename Ttype>
    void mergeAppendedData(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

    template<typename Ttype>
    Ttype* getAligned(size_t alignment, std::vector<Ttype> &data);

    template<typename I>
    constexpr I my_int_pow(I base, I exponent)
    {
        if(exponent < 0) return 0;
        if(exponent == 0) return 1;
        if(exponent == 1) return base;
        I tmp = my_int_pow<I>(base, exponent / 2);
        I result = tmp * tmp;
        if(exponent % 2 == 1) result *= base;
        return result;
    }

    template<typename T>
    constexpr size_t get_max_val_no_precision_loss_in_fp()
    {
        constexpr size_t result = my_int_pow<long long>(std::numeric_limits<T>::radix, std::numeric_limits<T>::digits);
        return result;
    }

    static void run_dummy_parallel_region();

    template<typename T> static constexpr bool is_real();
    template<> constexpr bool is_real<int32_t>() { return true; }
    template<> constexpr bool is_real<int64_t>() { return true; }
    template<> constexpr bool is_real<float>() { return true; }
    template<> constexpr bool is_real<double>() { return true; }
    template<> constexpr bool is_real<std::complex<float>>() { return false; }
    template<> constexpr bool is_real<std::complex<double>>() { return false; }

    template<typename T> static constexpr bool is_complex();
    template<> constexpr bool is_complex<int32_t>() { return false; }
    template<> constexpr bool is_complex<int64_t>() { return false; }
    template<> constexpr bool is_complex<float>() { return false; }
    template<> constexpr bool is_complex<double>() { return false; }
    template<> constexpr bool is_complex<std::complex<float>>() { return true; }
    template<> constexpr bool is_complex<std::complex<double>>() { return true; }

    template<typename T> static T & real_ref(std::complex<T> & x) { return reinterpret_cast<T*>(&x)[0]; }
    template<typename T> static T & imag_ref(std::complex<T> & x) { return reinterpret_cast<T*>(&x)[1]; }

    template<typename T> struct remove_complex { using type = decltype(std::real(T{})); };
    template<typename T> using remove_complex_t = typename remove_complex<T>::type;

    template<typename T>
    inline void atomic_add(T & dst, T val)
    {
        if constexpr(is_complex<T>()) {
            atomic_add(real_ref(dst), real_ref(val));
            atomic_add(imag_ref(dst), imag_ref(val));
        }
        else {
            #pragma omp atomic
            dst += val;
        }
    }

    template<typename I>
    I round_up(I num, I align)
    {
        static_assert(std::is_integral_v<I>);
        if constexpr(std::is_signed_v<I>) if(num == 0) return 0;
        return ((num - 1) / align + 1) * align;
    }
    
    template<typename I>
    I round_down(I num, I align)
    {
        static_assert(std::is_integral_v<I>);
        return (num / align) * align;
    }

    template<typename T>
    T * pointer_advance(T * ptr, size_t offset_B)
    {
        return (T*)((char*)ptr + offset_B);
    }

    template<typename ITER_IN, typename ITER_OUT, typename OP>
    inline void transform_n(ITER_IN begin, size_t n, ITER_OUT end, OP op)
    {
        std::transform(begin, begin + n, end, op);
    }
}

}

#include "utils.hpp"


#endif /* BASIS_UTILITIES_UTILS_H_ */
