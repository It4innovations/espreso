
#ifndef SRC_MATH_SIMD_SIMD_SSE2_H_
#define SRC_MATH_SIMD_SIMD_SSE2_H_

#if defined(__SSE2__)

#include "basis/utilities/inline.h"

#include <cstddef>
#include <cmath>
#include <immintrin.h>

#include <cstdio>

struct SIMD
{
	enum: size_t {
		size = 2U
	};

	ALWAYS_INLINE SIMD() noexcept: data(_mm_setzero_pd()) { }
	ALWAYS_INLINE SIMD(__m128d value) noexcept: data(value) { }
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
		__m128i tmp = _mm_set_epi32(1<<31, 0, 1<<31, 0);
		return _mm_xor_pd(data, reinterpret_cast<__m128d>(tmp));
	}

	ALWAYS_INLINE SIMD operator+ () const noexcept
	{
		return data;
	}

	__m128d data;
};

ALWAYS_INLINE const SIMD load1(const double &from) noexcept
{
	return _mm_load1_pd(&from);
}

ALWAYS_INLINE const SIMD load(const double *from) noexcept
{
	return _mm_load_pd(from);
}

ALWAYS_INLINE void store(double *to, const SIMD& value) noexcept
{
	_mm_store_pd(to, value.data);
}

ALWAYS_INLINE const SIMD operator+(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm_add_pd(v1.data, v2.data);
}

ALWAYS_INLINE const SIMD operator*(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm_mul_pd(v1.data, v2.data);
}

ALWAYS_INLINE const SIMD operator-(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm_sub_pd(v1.data, v2.data);
}

ALWAYS_INLINE const SIMD operator/(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm_div_pd(v1.data, v2.data);
}

ALWAYS_INLINE SIMD zeros() noexcept
{
	return _mm_setzero_pd();
}

ALWAYS_INLINE SIMD ones() noexcept
{
	return _mm_set1_pd(1.0);
}

ALWAYS_INLINE SIMD negate(const SIMD& value) noexcept
{
	__m128i tmp = _mm_set_epi32(1<<31, 0, 1<<31, 0);
	return _mm_xor_pd(value.data, reinterpret_cast<__m128d>(tmp));
}

ALWAYS_INLINE SIMD sqrt(const SIMD& value) noexcept
{
	return _mm_sqrt_pd(value.data);
}

ALWAYS_INLINE SIMD rsqrt14(const SIMD& v) noexcept // TODO: improve it
{
	__m128d sqrt = _mm_sqrt_pd(v.data);
	return __m128d{
		v.data[0] > 0. ? 1. / sqrt[0] : 0.,
		v.data[1] > 0. ? 1. / sqrt[1] : 0.
	};
}

ALWAYS_INLINE SIMD positive_guarded_recip(const SIMD& v) noexcept // TODO: improve it
{
	return __m128d{
		v.data[0] > 0. ? 1. / v.data[0] : 0.,
		v.data[1] > 0. ? 1. / v.data[1] : 0.
	};
}

ALWAYS_INLINE SIMD max(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm_max_pd(v1.data, v2.data);
}

ALWAYS_INLINE SIMD cos(const SIMD& value) noexcept
{
	return __m128d{
		std::cos(value.data[0]),
		std::cos(value.data[1])
	};
}

ALWAYS_INLINE SIMD acos(const SIMD& value) noexcept
{
	return __m128d{
		std::acos(value.data[0]),
		std::acos(value.data[1])
	};
}

ALWAYS_INLINE SIMD ispositive(const SIMD& v) noexcept
{
	return __m128d{
		v.data[0] > 0.0 ? 1.0 : 0.0,
		v.data[1] > 0.0 ? 1.0 : 0.0
	};
}

#endif // __SSE2__
#endif /* SRC_MATH_SIMD_SIMD_SSE2_H_ */
