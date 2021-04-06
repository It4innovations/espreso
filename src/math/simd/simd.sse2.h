
#ifndef SRC_MATH_SIMD_SIMD_SSE2_H_
#define SRC_MATH_SIMD_SIMD_SSE2_H_

#if defined(__SSE2__)

#include "basis/utilities/inline.h"

#include <cstddef>
#include <emmintrin.h>

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

	ALWAYS_INLINE double operator[](size_t i) const noexcept
	{
		return reinterpret_cast<const double*>(&data)[i];
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

ALWAYS_INLINE double sum(const SIMD& value) noexcept
{
	return _mm_cvtsd_f64(_mm_add_sd(value.data, _mm_unpackhi_pd(value.data, value.data)));
}

#endif // __SSE2__
#endif /* SRC_MATH_SIMD_SIMD_SSE2_H_ */
