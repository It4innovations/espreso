
#ifndef SRC_PHYSICS_ASSEMBLER_MATH_HPP_
#define SRC_PHYSICS_ASSEMBLER_MATH_HPP_

namespace espreso {

// |cos , -sin| |c[0] , c[2]| | cos , sin|
// |sin ,  cos| |c[3] , c[1]| |-sin , cos|
// |cos * c[0] - sin * c[3] , cos * c[2] - sin * c[1]|  | cos , sin|
// |sin * c[0] + cos * c[3] , sin * c[2] + cos * c[1]|  |-sin , cos|
static inline void rotate2x2(const double &cos, const double &sin, double * __restrict__ m22)
{
	double origin[4] = { m22[0], m22[1], m22[2], m22[3] };
	m22[0] = (cos * origin[0] - sin * origin[2]) * cos - (cos * origin[1] - sin * origin[3]) * sin;
	m22[1] = (cos * origin[0] - sin * origin[2]) * sin + (cos * origin[1] - sin * origin[3]) * cos;
	m22[2] = (sin * origin[0] + cos * origin[2]) * cos - (sin * origin[1] + cos * origin[3]) * sin;
	m22[3] = (sin * origin[0] + cos * origin[2]) * sin + (sin * origin[1] + cos * origin[3]) * cos;
}

static inline void rotate3x3(const double * __restrict__ cos, const double * __restrict__ sin, double * __restrict__ m33)
{
	double t[3][3] {
		{ cos[1] * cos[2]                           , cos[1] * sin[2]                           , -sin[1] },
		{ cos[2] * sin[0] * sin[1] - cos[0] * sin[2], cos[0] * cos[2] + sin[0] * sin[1] * sin[2], cos[1] * sin[0] },
		{ sin[0] * sin[2] + cos[0] * cos[2] * sin[1], cos[0] * sin[1] * sin[2] - cos[2] * sin[0], cos[0] * cos[1] }
	};

	double origin[9] = { m33[0], m33[1], m33[2], m33[3], m33[4], m33[5], m33[6], m33[7], m33[8] };

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			double _a = t[0][i] * origin[0] + t[1][i] * origin[3] + t[2][i] * origin[6];
			double _b = t[0][i] * origin[1] + t[1][i] * origin[4] + t[2][i] * origin[7];
			double _c = t[0][i] * origin[2] + t[1][i] * origin[5] + t[2][i] * origin[8];
			m33[3 * i + j] = _a * t[0][j] + _b * t[1][j] + _c * t[2][j];
		}
	}
}

template<int nodes, int dimension>
static void NtoGP(const double * __restrict__ gp_N, const double * __restrict__ nvalues, double * __restrict__ gpvalues)
{
	for (int d = 0; d < dimension; ++d) {
		gpvalues[d] = 0;
	}
	for (int n = 0; n < nodes; ++n) {
		for (int d = 0; d < dimension; ++d) {
			gpvalues[d] += gp_N[n] * nvalues[dimension * n + d];
		}
	}
}

template<int nodes, int dimension>
static void NtoGPSimd(const double * __restrict__ gp_N, const double * __restrict__ nvalues, double * __restrict__ gpvalues)
{

	SIMD res[dimension];

	for (int d = 0; d < dimension; ++d) {
		res[d] = zeros();
	}

	for (int n = 0; n < nodes; ++n) {
		SIMD gp = load(&gp_N[n * SIMD::size]);

		for (int d = 0; d < dimension; ++d) {
			SIMD nv = load(&nvalues[(dimension * n + d) * SIMD::size]);
			res[d] = res[d] + gp * nv;
		}
	}

	for (int d = 0; d < dimension; ++d) {
		store(&gpvalues[d * SIMD::size], res[d]);
	}
}

template<int N>
static inline void M12M2N(const double * __restrict__ m12, const double * __restrict__ m2N, double * __restrict__ result)
{
	for (int n = 0; n < N; ++n) {
		result[0 + n] = m12[0] * m2N[n] + m12[1] * m2N[N + n];
	}
}

template<int N>
static inline void M13M3N(const double * __restrict__ m13, const double * __restrict__ m3N, double * __restrict__ result)
{
	for (int n = 0; n < N; ++n) {
		result[n] = m13[0] * m3N[n] + m13[1] * m3N[N + n] + m13[2] * m3N[2 * N + n];
	}
}

template<int N>
static inline void M22M2N(const double * __restrict__ m22, const double * __restrict__ m2N, double * __restrict__ result)
{
	for (int n = 0; n < N; ++n) {
		result[0 + n] = m22[0] * m2N[n] + m22[1] * m2N[N + n];
		result[N + n] = m22[2] * m2N[n] + m22[3] * m2N[N + n];
	}
}

template<int N>
static inline void M33M3N(const double * __restrict__ m33, const double * __restrict__ m3N, double * __restrict__ result)
{
	for (int n = 0; n < N; ++n) {
		result[0 * N + n] = m33[0] * m3N[n] + m33[1] * m3N[N + n] + m33[2] * m3N[2 * N + n];
		result[1 * N + n] = m33[3] * m3N[n] + m33[4] * m3N[N + n] + m33[5] * m3N[2 * N + n];
		result[2 * N + n] = m33[6] * m3N[n] + m33[7] * m3N[N + n] + m33[8] * m3N[2 * N + n];
	}
}

template<int N>
static inline void MN1M1N(const double &sumscale, const double * __restrict__ mN1, const double * __restrict__ m1N, double * __restrict__ mNN)
{
	for (int n = 0; n < N; ++n) {
		for (int m = 0; m < N; ++m) {
			mNN[n * N + m] = sumscale * mN1[n] * m1N[m];
		}
	}
}

template<int N>
static inline void M1NMN2(const double &sumscale, const double * __restrict__ m1N, const double * __restrict__ mN2, double * __restrict__ m12)
{
	m12[2 * 0 + 0] = 0;
	m12[2 * 0 + 1] = 0;
	for (int n = 0; n < N; ++n) {
		m12[2 * 0 + 0] += m1N[N * 0 + n] * mN2[2 * n + 0];
		m12[2 * 0 + 1] += m1N[N * 0 + n] * mN2[2 * n + 1];
	}
	m12[2 * 0 + 0] *= sumscale;
	m12[2 * 0 + 1] *= sumscale;
}

template<int N>
static inline void M1NMN3(const double &sumscale, const double * __restrict__ m1N, const double * __restrict__ mN3, double * __restrict__ m13)
{
	m13[3 * 0 + 0] = 0;
	m13[3 * 0 + 1] = 0;
	m13[3 * 0 + 2] = 0;
	for (int n = 0; n < N; ++n) {
		m13[3 * 0 + 0] += m1N[N * 0 + n] * mN3[3 * n + 0];
		m13[3 * 0 + 1] += m1N[N * 0 + n] * mN3[3 * n + 1];
		m13[3 * 0 + 2] += m1N[N * 0 + n] * mN3[3 * n + 2];
	}
	m13[3 * 0 + 0] *= sumscale;
	m13[3 * 0 + 1] *= sumscale;
	m13[3 * 0 + 2] *= sumscale;
}

template<int N>
static inline void M2NMN3(const double &sumscale, const double * __restrict__ m2N, const double * __restrict__ mN3, double * __restrict__ m23)
{
	m23[3 * 0 + 0] = 0;
	m23[3 * 0 + 1] = 0;
	m23[3 * 0 + 2] = 0;
	m23[3 * 1 + 0] = 0;
	m23[3 * 1 + 1] = 0;
	m23[3 * 1 + 2] = 0;
	for (int n = 0; n < N; ++n) {
		m23[3 * 0 + 0] += m2N[N * 0 + n] * mN3[3 * n + 0];
		m23[3 * 0 + 1] += m2N[N * 0 + n] * mN3[3 * n + 1];
		m23[3 * 0 + 2] += m2N[N * 0 + n] * mN3[3 * n + 2];
		m23[3 * 1 + 0] += m2N[N * 1 + n] * mN3[3 * n + 0];
		m23[3 * 1 + 1] += m2N[N * 1 + n] * mN3[3 * n + 1];
		m23[3 * 1 + 2] += m2N[N * 1 + n] * mN3[3 * n + 2];
	}
	m23[3 * 0 + 0] *= sumscale;
	m23[3 * 0 + 1] *= sumscale;
	m23[3 * 0 + 2] *= sumscale;
	m23[3 * 1 + 0] *= sumscale;
	m23[3 * 1 + 1] *= sumscale;
	m23[3 * 1 + 2] *= sumscale;
}

template<int N>
static inline void ADDM2NMN1(const double &sumscale, const double * __restrict__ m2N, const double * __restrict__ mN1, double * __restrict__ m21)
{
	double x[2] = { 0, 0 };
	for (int n = 0; n < N; ++n) {
		x[0] += m2N[0 * N + n] * mN1[n];
		x[1] += m2N[1 * N + n] * mN1[n];
	}
	m21[0] += sumscale * x[0];
	m21[1] += sumscale * x[1];
}

template<int N>
static inline void ADDM3NMN1(const double &sumscale, const double * __restrict__ m3N, const double * __restrict__ mN1, double * __restrict__ m31)
{
	double x[3] = { 0, 0, 0 };
	for (int n = 0; n < N; ++n) {
		x[0] += m3N[0 * N + n] * mN1[n];
		x[1] += m3N[1 * N + n] * mN1[n];
		x[2] += m3N[2 * N + n] * mN1[n];
	}
	m31[0] += sumscale * x[0];
	m31[1] += sumscale * x[1];
	m31[2] += sumscale * x[2];
}

template<int N>
static inline void ADDM22M2NMN1(const double &sumscale, const double * __restrict__ m22, const double * __restrict__ m2N, const double * __restrict__ mN1, double * __restrict__ m21)
{
	double x[2] = { 0, 0 };
	for (int n = 0; n < N; ++n) {
		x[0] += m2N[0 * N + n] * mN1[n];
		x[1] += m2N[1 * N + n] * mN1[n];
	}
	m21[0] += sumscale * (m22[0 * 2 + 0] * x[0] + m22[0 * 2 + 1] * x[1]);
	m21[1] += sumscale * (m22[1 * 2 + 0] * x[0] + m22[1 * 2 + 1] * x[1]);
}

template<int N>
static inline void ADDM33M3NMN1(const double &sumscale, const double * __restrict__ m33, const double * __restrict__ m3N, const double * __restrict__ mN1, double * __restrict__ m31)
{
	double x[3] = { 0, 0, 0 };
	for (int n = 0; n < N; ++n) {
		x[0] += m3N[0 * N + n] * mN1[n];
		x[1] += m3N[1 * N + n] * mN1[n];
		x[2] += m3N[2 * N + n] * mN1[n];
	}
	m31[0] += sumscale * (m33[0 * 3 + 0] * x[0] + m33[0 * 3 + 1] * x[1] + m33[0 * 3 + 2] * x[2]);
	m31[1] += sumscale * (m33[1 * 3 + 0] * x[0] + m33[1 * 3 + 1] * x[1] + m33[1 * 3 + 2] * x[2]);
	m31[2] += sumscale * (m33[2 * 3 + 0] * x[0] + m33[2 * 3 + 1] * x[1] + m33[2 * 3 + 2] * x[2]);
}

template<int N>
static inline void ADDMN1M1N(const double &sumscale, const double * __restrict__ mN1, const double * __restrict__ m1N, double * __restrict__ mNN)
{
	for (int n = 0; n < N; ++n) {
		for (int m = 0; m < N; ++m) {
			mNN[n * N + m] += sumscale * (mN1[n] * m1N[m]);
		}
	}
}

template<int N>
static inline void ADDMN1M1N(const double &sumscale, const double * __restrict__ m1N, double * __restrict__ mNN)
{
	for (int n = 0; n < N; ++n) {
		for (int m = 0; m < N; ++m) {
			mNN[n * N + m] += sumscale * (m1N[n] * m1N[m]);
		}
	}
}

template<int N>
static inline void ADDMN2M2N(const double &sumscale, const double * __restrict__ m2N, double * __restrict__ mNN)
{
	for (int n = 0; n < N; ++n) {
		for (int m = 0; m < N; ++m) {
			mNN[n * N + m] += sumscale * (m2N[n] * m2N[m] + m2N[N + n] * m2N[N + m]);
		}
	}
}

template<int N>
static inline void ADDMN3M3N(const double &sumscale, const double * __restrict__ m3N, double * __restrict__ mNN)
{
	for (int n = 0; n < N; ++n) {
		for (int m = 0; m < N; ++m) {
			mNN[n * N + m] += sumscale * (m3N[n] * m3N[m] + m3N[N + n] * m3N[N + m] + m3N[2 * N + n] * m3N[2 * N + m]);
		}
	}
}

template<int N>
static inline void ADDMN2M22M2N(const double &sumscale, const double * __restrict__ m22, const double * __restrict__ m2N, double * __restrict__ mNN)
{
	for (int n = 0; n < N; ++n) {
		for (int m = 0; m < N; ++m) {
			double a = m2N[n] * m22[0] + m2N[N + n] * m22[2];
			double b = m2N[n] * m22[1] + m2N[N + n] * m22[3];
			mNN[n * N + m] += sumscale * (a * m2N[m] + b * m2N[N + m]);
		}
	}
}

template<int N>
static inline void ADDMN3M33M3N(const double &sumscale, const double * __restrict__ m33, const double * __restrict__ m3N, double * __restrict__ mNN)
{
	for (int n = 0; n < N; ++n) {
		for (int m = 0; m < N; ++m) {
			double a = m3N[n] * m33[0] + m3N[N + n] * m33[3] + m3N[2 * N + n] * m33[6];
			double b = m3N[n] * m33[1] + m3N[N + n] * m33[4] + m3N[2 * N + n] * m33[7];
			double c = m3N[n] * m33[2] + m3N[N + n] * m33[5] + m3N[2 * N + n] * m33[8];
			mNN[n * N + m] += sumscale * (a * m3N[m] + b * m3N[N + m] + c * m3N[2 * N + m]);
		}
	}
}

}

#endif /* SRC_PHYSICS_ASSEMBLER_MATH_HPP_ */
