
#ifndef SRC_MATH_SIMD_SIMD_H_
#define SRC_MATH_SIMD_SIMD_H_

#if defined(__AVX__)
#include "simd.avx.h"
#elif defined(__SSE2__)
#include "simd.sse2.h"
#else

#include "basis/utilities/inline.h"

#include <cstddef>

struct SIMD
{
	enum: size_t {
		size = 1U
	};

	ALWAYS_INLINE SIMD() noexcept: data(0) { }
	ALWAYS_INLINE SIMD(double value) noexcept: data(value) { }
	ALWAYS_INLINE SIMD(const SIMD &other) noexcept: data(other.data) { }

	ALWAYS_INLINE SIMD& operator=(const SIMD &other) noexcept
	{
		data = other.data;
		return *this;
	}

	ALWAYS_INLINE double operator[](size_t i) const noexcept
	{
		return data;
	}

	double data;
};

ALWAYS_INLINE const SIMD load(const double *from) noexcept
{
	return *from;
}

ALWAYS_INLINE void store(double *to, const SIMD& value) noexcept
{
	*to = value.data;
}

ALWAYS_INLINE const SIMD operator+(const SIMD& v1, const SIMD& v2) noexcept
{
	return v1.data +  v2.data;
}

ALWAYS_INLINE const SIMD operator*(const SIMD& v1, const SIMD& v2) noexcept
{
	return v1.data *  v2.data;
}

ALWAYS_INLINE double sum(const SIMD& value) noexcept
{
	return value.data;
}

#endif // default SIMD with single double
#endif /* SRC_MATH_SIMD_SIMD_H_ */
