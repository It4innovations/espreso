
#ifndef SRC_MATH_SIMD_SIMD_AVX_H_
#define SRC_MATH_SIMD_SIMD_AVX_H_

#if defined(__AVX__)

#include "basis/utilities/inline.h"

#include <cstddef>
#include <cmath>
#include <immintrin.h>

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
		__m256i tmp = _mm256_set_epi32(1<<31, 0, 1<<31, 0, 1<<31, 0, 1<<31, 0);
		return _mm256_xor_pd(data, reinterpret_cast<__m256d>(tmp));
	}

	ALWAYS_INLINE SIMD operator+ () const noexcept
	{
		return data;
	}

	__m256d data;
};

ALWAYS_INLINE const SIMD load1(const double &from) noexcept
{
	return _mm256_broadcast_sd(&from);
}

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

ALWAYS_INLINE const SIMD operator-(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm256_sub_pd(v1.data, v2.data);
}

ALWAYS_INLINE const SIMD operator/(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm256_div_pd(v1.data, v2.data);
}

ALWAYS_INLINE SIMD zeros() noexcept
{
	return _mm256_setzero_pd();
}

ALWAYS_INLINE SIMD ones() noexcept
{
	return _mm256_set1_pd(1.0);
}
ALWAYS_INLINE SIMD negate(const SIMD& value) noexcept
{
	__m256i tmp = _mm256_set_epi32(1<<31, 0, 1<<31, 0, 1<<31, 0, 1<<31, 0);
	return _mm256_xor_pd(value.data, reinterpret_cast<__m256d>(tmp));
}

ALWAYS_INLINE SIMD sqrt(const SIMD& value) noexcept
{
	return _mm256_sqrt_pd(value.data);
}

ALWAYS_INLINE SIMD rsqrt14(const SIMD& v) noexcept
{
	__m256d sqrt = _mm256_sqrt_pd(v.data);
	return __m256d{
		v.data[0] > 0. ? 1. / sqrt[0] : 0.,
		v.data[1] > 0. ? 1. / sqrt[1] : 0.,
		v.data[2] > 0. ? 1. / sqrt[2] : 0.,
		v.data[3] > 0. ? 1. / sqrt[3] : 0.
	};
}

ALWAYS_INLINE SIMD log(const SIMD& v) noexcept
{
    return __m256d{
        std::log(v.data[0]),
        std::log(v.data[1]),
        std::log(v.data[2]),
        std::log(v.data[3])
    };
}

ALWAYS_INLINE SIMD exp(const SIMD& v) noexcept
{
    return __m256d{
        std::exp(v.data[0]),
        std::exp(v.data[1]),
        std::exp(v.data[2]),
        std::exp(v.data[3])
    };
}

ALWAYS_INLINE SIMD abs(const SIMD& v) noexcept
{
    return __m256d{
        std::abs(v.data[0]),
        std::abs(v.data[1]),
        std::abs(v.data[2]),
        std::abs(v.data[3])
    };
}

ALWAYS_INLINE SIMD positive_guarded_recip(const SIMD& v) noexcept // TODO: improve it
{
	return __m256d{
		v.data[0] > 0. ? 1. / v.data[0] : 0.,
		v.data[1] > 0. ? 1. / v.data[1] : 0.,
		v.data[2] > 0. ? 1. / v.data[2] : 0.,
		v.data[3] > 0. ? 1. / v.data[3] : 0.
	};
}

ALWAYS_INLINE SIMD max(const SIMD& v1, const SIMD& v2) noexcept
{
	return _mm256_max_pd(v1.data, v2.data);
}

ALWAYS_INLINE SIMD cos(const SIMD& value) noexcept
{
	return __m256d{
		std::cos(value.data[0]),
		std::cos(value.data[1]),
		std::cos(value.data[2]),
		std::cos(value.data[3])
	};
}

ALWAYS_INLINE SIMD acos(const SIMD& value) noexcept
{
	return __m256d{
		std::acos(value.data[0]),
		std::acos(value.data[1]),
		std::acos(value.data[2]),
		std::acos(value.data[3])
	};
}

ALWAYS_INLINE SIMD ispositive(const SIMD& v) noexcept
{
	return __m256d{
		v.data[0] > 0.0 ? 1.0 : 0.0,
		v.data[1] > 0.0 ? 1.0 : 0.0,
		v.data[2] > 0.0 ? 1.0 : 0.0,
		v.data[3] > 0.0 ? 1.0 : 0.0
	};
}

#endif // __AVX__
#endif /* SRC_MATH_SIMD_SIMD_AVX_H_ */
