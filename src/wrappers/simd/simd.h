
#ifndef SRC_MATH_SIMD_SIMD_H_
#define SRC_MATH_SIMD_SIMD_H_

#if !defined(SIMD_OFF) && defined(__AVX512F__) && defined(__AVX512DQ__)
#include "simd.avx512.h"
#elif !defined(SIMD_OFF) && defined(__AVX__)
#include "simd.avx.h"
#elif !defined(SIMD_OFF) && defined(__SSE2__)
#include "simd.sse2.h"
#elif !defined(SIMD_OFF) && defined(__ARM_FEATURE_SVE) && defined(SIMD_ARM_SVE_GENERAL)
#include "simd.sve.h"
#elif !defined(SIMD_OFF) && defined(__ARM_FEATURE_SVE) && defined(SIMD_ARM_SVE_ARRAY512)
#include "simd.sve.array512.h"
#else

#include "basis/utilities/inline.h"

#include <cstddef>
#include <cmath>

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

	ALWAYS_INLINE double& operator[](size_t i) noexcept
	{
		return reinterpret_cast<double*>(&data)[i];
	}

	ALWAYS_INLINE const double& operator[](size_t i) const noexcept
	{
		return reinterpret_cast<const double*>(&data)[i];
	}

	ALWAYS_INLINE double operator- () const noexcept
	{
		return -data;
	}

	ALWAYS_INLINE double operator+ () const noexcept
	{
		return data;
	}

	double data;
};

ALWAYS_INLINE const SIMD load1(const double &from) noexcept
{
	return from;
}

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
	return v1.data + v2.data;
}

ALWAYS_INLINE const SIMD operator*(const SIMD& v1, const SIMD& v2) noexcept
{
	return v1.data * v2.data;
}

ALWAYS_INLINE const SIMD operator-(const SIMD& v1, const SIMD& v2) noexcept
{
	return v1.data - v2.data;
}

ALWAYS_INLINE const SIMD operator/(const SIMD& v1, const SIMD& v2) noexcept
{
	return v1.data / v2.data;
}

ALWAYS_INLINE SIMD zeros() noexcept
{
	return 0.0;
}

ALWAYS_INLINE SIMD ones() noexcept
{
	return 1.0;
}

ALWAYS_INLINE SIMD negate(const SIMD& value) noexcept
{
	return -value.data;
}

ALWAYS_INLINE SIMD sqrt(const SIMD& value) noexcept
{
	return std::sqrt(value.data);
}

ALWAYS_INLINE SIMD rsqrt14(const SIMD& v) noexcept
{
	if (v.data > .0) {
		return 1. / std::sqrt(v.data);
	}
	return 0;
}

ALWAYS_INLINE SIMD log(const SIMD& value) noexcept
{
    return std::log(value.data);
}

ALWAYS_INLINE SIMD positive_guarded_recip(const SIMD& v) noexcept // TODO: improve it
{
	return v.data > 0. ? 1. / v.data : 0.;
}

ALWAYS_INLINE SIMD max(const SIMD& v1, const SIMD& v2) noexcept
{
	return std::max(v1.data, v2.data);
}

ALWAYS_INLINE SIMD cos(const SIMD& v) noexcept
{
	return std::cos(v.data);
}

ALWAYS_INLINE SIMD acos(const SIMD& v) noexcept
{
	return std::acos(v.data);
}

ALWAYS_INLINE SIMD ispositive(const SIMD& v) noexcept
{
	return v.data > 0.0 ? 1.0 : 0.0;
}
#endif // default SIMD with single double
#endif /* SRC_MATH_SIMD_SIMD_H_ */
