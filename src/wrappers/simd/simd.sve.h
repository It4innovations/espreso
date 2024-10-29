
#ifndef SRC_MATH_SIMD_SIMD_SVE_H_
#define SRC_MATH_SIMD_SIMD_SVE_H_

#if defined(__ARM_FEATURE_SVE)

#include "basis/utilities/inline.h"

#include <cstddef>
#include <cmath>
#include <arm_sve.h>

#if __ARM_FEATURE_SVE_BITS==0
#error "Set __ARM_FEATURE_SVE_BITS"
#endif

typedef svfloat64_t __sved __attribute__((arm_sve_vector_bits(512)));
typedef svbool_t    __svep __attribute__((arm_sve_vector_bits(512)));

struct SIMD
{
	enum: size_t {
		size = __ARM_FEATURE_SVE_BITS / 64
	};

	ALWAYS_INLINE SIMD() noexcept: data(svdup_n_f64(0.0)) { }
	ALWAYS_INLINE SIMD(__sved value) noexcept: data(value) { }
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
		return svneg_f64_x(svptrue_b64(), data);
	}

	ALWAYS_INLINE SIMD operator+ () const noexcept
	{
		return data;
	}

	__sved data;
};

ALWAYS_INLINE const SIMD load1(const double &from) noexcept
{
	return svdup_n_f64(from);
}

ALWAYS_INLINE const SIMD load(const double *from) noexcept
{
	return svld1_f64(svptrue_b64(), from);
}

ALWAYS_INLINE void store(double *to, const SIMD& value) noexcept
{
	svst1_f64(svptrue_b64(), to, value.data);
}

ALWAYS_INLINE const SIMD operator+(const SIMD& v1, const SIMD& v2) noexcept
{
	return svadd_f64_x(svptrue_b64(), v1.data, v2.data);
}

ALWAYS_INLINE const SIMD operator*(const SIMD& v1, const SIMD& v2) noexcept
{
	return svmul_f64_x(svptrue_b64(), v1.data, v2.data);
}

ALWAYS_INLINE const SIMD operator-(const SIMD& v1, const SIMD& v2) noexcept
{
	return svsub_f64_x(svptrue_b64(), v1.data, v2.data);
}

ALWAYS_INLINE const SIMD operator/(const SIMD& v1, const SIMD& v2) noexcept
{
	return svdiv_f64_x(svptrue_b64(), v1.data, v2.data);
}

ALWAYS_INLINE SIMD zeros() noexcept
{
	return svdup_n_f64(0.0);
}

ALWAYS_INLINE SIMD ones() noexcept
{
	return svdup_n_f64(1.0);
}

ALWAYS_INLINE SIMD negate(const SIMD& value) noexcept
{
	return svneg_f64_x(svptrue_b64(), value.data);
}

ALWAYS_INLINE SIMD sqrt(const SIMD& value) noexcept
{
	return svsqrt_f64_x(svptrue_b64(), value.data);
}

ALWAYS_INLINE SIMD rsqrt14(const SIMD& value) noexcept // TODO: improve it
{
	return svrsqrte_f64(value.data);
}

ALWAYS_INLINE SIMD log(const SIMD& v) noexcept
{
    return __sved{};
}

ALWAYS_INLINE SIMD exp(const SIMD& v) noexcept
{
    return __sved{};
}

ALWAYS_INLINE SIMD abs(const SIMD& v) noexcept
{
    return __sved{};
}


ALWAYS_INLINE SIMD positive_guarded_recip(const SIMD& value) noexcept // TODO: improve it
{
	return svrsqrte_f64(value.data);
//	return __sved{
//		v.data[0] > 0. ? 1. / v.data[0] : 0.,
//		v.data[1] > 0. ? 1. / v.data[1] : 0.
//	};
}

ALWAYS_INLINE SIMD max(const SIMD& v1, const SIMD& v2) noexcept
{
	return svmax_f64_x(svptrue_b64(), v1.data, v2.data);
}

ALWAYS_INLINE SIMD cos(const SIMD& value) noexcept
{
	return __sved{};
//	return __sved{
//		std::cos(value.data[0]),
//		std::cos(value.data[1])
//	};
}

ALWAYS_INLINE SIMD acos(const SIMD& value) noexcept
{
	return __sved{};
//	return __sved{
//		std::acos(value.data[0]),
//		std::acos(value.data[1])
//	};
}

ALWAYS_INLINE SIMD ispositive(const SIMD& v) noexcept
{
	return __sved{};
//	return __sved{
//		v.data[0] > 0.0 ? 1.0 : 0.0,
//		v.data[1] > 0.0 ? 1.0 : 0.0
//	};
}



#endif /* __ARM_FEATURE_SVE */
#endif /* SRC_MATH_SIMD_SIMD_SVE_H_ */
