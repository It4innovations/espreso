
#ifndef SRC_MATH_SIMD_SIMD_SSE2_H_
#define SRC_MATH_SIMD_SIMD_SSE2_H_

#if defined(__SSE2__)

#include "basis/utilities/inline.h"

#include <cstddef>
#include <immintrin.h>

struct SIMD
{
	enum: size_t {
		size = 2U
	};

	ALWAYS_INLINE SIMD() noexcept: data(_mm_setzero_pd()) { }
	ALWAYS_INLINE SIMD(__m128d value) noexcept: data(value) { }
	ALWAYS_INLINE SIMD(const SIMD &other) noexcept: data(other.data) { }

	ALWAYS_INLINE void fill(const double &value) noexcept
	{
		data[0] = data[1] = value;
	}

	ALWAYS_INLINE operator double() const noexcept
	{
		return data[0];
	}

	ALWAYS_INLINE operator double*() noexcept
	{
		return reinterpret_cast<double*>(&data);
	}

	ALWAYS_INLINE SIMD& operator=(const double &value) noexcept
	{
		data[0] = value;
		return *this;
	}

	ALWAYS_INLINE SIMD& operator=(const SIMD &other) noexcept
	{
		data = other.data;
		return *this;
	}

	ALWAYS_INLINE double& operator[](size_t i) noexcept
	{
		return reinterpret_cast<double*>(&data)[i];
	}


	ALWAYS_INLINE double operator[](size_t i) const noexcept
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

#endif // __SSE2__
#endif /* SRC_MATH_SIMD_SIMD_SSE2_H_ */
