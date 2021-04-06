
#ifndef SRC_MATH_SIMD_SIMD_AVX_H_
#define SRC_MATH_SIMD_SIMD_AVX_H_

#if defined(__AVX__)

#include "basis/utilities/inline.h"

#include <cstddef>
#include <emmintrin.h>

struct SIMD
{
	enum: size_t {
		size = 4U
	};

	ALWAYS_INLINE SIMD() noexcept: data(_mm256_setzero_pd()) { }
	ALWAYS_INLINE SIMD(__m256d value) noexcept: data(value) { }
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

	__m256d data;
};

ALWAYS_INLINE const SIMD load(const double *from) noexcept
{
	return _mm256_load_pd(from);
}

ALWAYS_INLINE void store(double *to, const SIMD& value) noexcept
{
	_mm256_store_pd(to, value.data);
}

ALWAYS_INLINE const SIMD operator+(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm256_add_pd(v1.data, v2.data);
}

ALWAYS_INLINE const SIMD operator*(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm256_mul_pd(v1.data, v2.data);
}

ALWAYS_INLINE double sum(const SIMD& value) noexcept
{
	const __m128d v128(_mm_add_pd(_mm256_castpd256_pd128(value.value), _mm256_extractf128_pd(value.value, 1)));
	return _mm_cvtsd_f64(_mm_add_sd(v128, _mm_unpackhi_pd(v128, v128)));
}

#endif // __AVX__
#endif /* SRC_MATH_SIMD_SIMD_AVX_H_ */
