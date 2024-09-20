
#ifndef SRC_MATH_SIMD_SIMD_AVX_H_
#define SRC_MATH_SIMD_SIMD_AVX_H_

#if defined(__AVX512F__) && defined(__AVX512DQ__)

#include "basis/utilities/inline.h"

#include <cstddef>
#include <cmath>
#include <immintrin.h>

struct SIMD
{
	enum: size_t {
		size = 8U
	};

	ALWAYS_INLINE SIMD() noexcept: data(_mm512_setzero_pd()) { }
	ALWAYS_INLINE SIMD(__m512d value) noexcept: data(value) { }
	ALWAYS_INLINE SIMD(const SIMD &other) noexcept: data(other.data) { }

	ALWAYS_INLINE SIMD& operator=(const SIMD &other) noexcept
	{
		data = other.data;
		return *this;
	}

	ALWAYS_INLINE double& operator[](size_t i) noexcept
	{
		return reinterpret_cast<double*>(&data)[i];
	}

	ALWAYS_INLINE const double& operator[](size_t i) const noexcept
	{
		return reinterpret_cast<const double*>(&data)[i];
	}
	
	ALWAYS_INLINE SIMD operator- () const noexcept
	{
		__m512i tmp = _mm512_set_epi32(1<<31, 0, 1<<31, 0, 1<<31, 0, 1<<31, 0, 1<<31, 0, 1<<31, 0, 1<<31, 0, 1<<31, 0);
		return _mm512_xor_pd(data, reinterpret_cast<__m512d>(tmp));
	}
	
	ALWAYS_INLINE SIMD operator+ () const noexcept
	{
		return data;
	}

	__m512d data;
};

ALWAYS_INLINE const SIMD load1(const double &from) noexcept
{
	return _mm512_set1_pd(from);
}

ALWAYS_INLINE const SIMD load(const double *from) noexcept
{
	return _mm512_load_pd(from);
}

ALWAYS_INLINE void store(double *to, const SIMD& value) noexcept
{
	_mm512_store_pd(to, value.data);
}

ALWAYS_INLINE const SIMD operator+(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm512_add_pd(v1.data, v2.data);
}

ALWAYS_INLINE const SIMD operator*(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm512_mul_pd(v1.data, v2.data);
}

ALWAYS_INLINE const SIMD operator-(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm512_sub_pd(v1.data, v2.data);
}

ALWAYS_INLINE const SIMD operator/(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm512_div_pd(v1.data, v2.data);
}

ALWAYS_INLINE SIMD zeros() noexcept
{
	return _mm512_setzero_pd();
}

ALWAYS_INLINE SIMD ones() noexcept
{
	return _mm512_set1_pd(1.0);
}
ALWAYS_INLINE SIMD negate(const SIMD& value) noexcept
{
	__m512i tmp = _mm512_set_epi32(1<<31, 0, 1<<31, 0, 1<<31, 0, 1<<31, 0, 1<<31, 0, 1<<31, 0, 1<<31, 0, 1<<31, 0);
	return _mm512_xor_pd(value.data, reinterpret_cast<__m512d>(tmp));
}

ALWAYS_INLINE SIMD sqrt(const SIMD& value) noexcept
{
	return _mm512_sqrt_pd(value.data);
}

ALWAYS_INLINE SIMD rsqrt14(const SIMD& v) noexcept
{
	__m512d zero = { 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14 };
	__mmask8 mask = _mm512_cmp_pd_mask(zero, v.data, _CMP_NLT_UQ);
	return _mm512_maskz_rsqrt14_pd(mask, v.data);
}

ALWAYS_INLINE SIMD log(const SIMD& v) noexcept
{
    return __m512d{
        std::log(v.data[0]),
        std::log(v.data[1]),
        std::log(v.data[2]),
        std::log(v.data[3]),
        std::log(v.data[4]),,
        std::log(v.data[5]),
        std::log(v.data[6]),
        std::log(v.data[7])
    };
}

ALWAYS_INLINE SIMD positive_guarded_recip(const SIMD& v) noexcept // TODO: improve it
{
	return __m512d{
		v.data[0] > 0. ? 1. / v.data[0] : 0.,
		v.data[1] > 0. ? 1. / v.data[1] : 0.,
		v.data[2] > 0. ? 1. / v.data[2] : 0.,
		v.data[3] > 0. ? 1. / v.data[3] : 0.,
		v.data[4] > 0. ? 1. / v.data[4] : 0.,
		v.data[5] > 0. ? 1. / v.data[5] : 0.,
		v.data[6] > 0. ? 1. / v.data[6] : 0.,
		v.data[7] > 0. ? 1. / v.data[7] : 0.
	};
}

ALWAYS_INLINE SIMD max(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm512_max_pd(v1.data, v2.data);
}

ALWAYS_INLINE SIMD cos(const SIMD& value) noexcept
{
	return __m512d{
			std::cos(value.data[0]),
			std::cos(value.data[1]),
			std::cos(value.data[2]),
			std::cos(value.data[3]),
			std::cos(value.data[4]),
			std::cos(value.data[5]),
			std::cos(value.data[6]),
			std::cos(value.data[7])
		};
}

ALWAYS_INLINE SIMD acos(const SIMD& value) noexcept
{
	return __m512d{
			std::acos(value.data[0]),
			std::acos(value.data[1]),
			std::acos(value.data[2]),
			std::acos(value.data[3]),
			std::acos(value.data[4]),
			std::acos(value.data[5]),
			std::acos(value.data[6]),
			std::acos(value.data[7])
		};
}

ALWAYS_INLINE SIMD ispositive(const SIMD& v) noexcept
{
	return __m512d{
		v.data[0] > 0.0 ? 1.0 : 0.0,
		v.data[1] > 0.0 ? 1.0 : 0.0,
		v.data[2] > 0.0 ? 1.0 : 0.0,
		v.data[3] > 0.0 ? 1.0 : 0.0,
		v.data[4] > 0.0 ? 1.0 : 0.0,
		v.data[5] > 0.0 ? 1.0 : 0.0,
		v.data[6] > 0.0 ? 1.0 : 0.0,
		v.data[7] > 0.0 ? 1.0 : 0.0
	};
}

#endif // __AVX__
#endif /* SRC_MATH_SIMD_SIMD_AVX_H_ */
