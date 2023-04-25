
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_GRADIENT_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_GRADIENT_H_

#include "subkernel.h"

namespace espreso {

struct TemperatureGradient: SubKernel {
	const char* name() const { return "TemperatureGradient"; }

	double *gradient, *end;

	TemperatureGradient()
	: gradient(nullptr), end(nullptr)
	{
		isconst = false;
		action = Assembler::SOLUTION;
	}

	void activate(size_t interval, NamedData *gradient)
	{
		this->gradient = gradient->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
		this->end = gradient->data.data() + gradient->data.size();
		isactive = 1;
	}
};

template <size_t nodes, size_t gps, size_t ndim, class Physics> struct TemperatureGradientKernel;

template <size_t nodes, size_t gps, class Physics>
struct TemperatureGradientKernel<nodes, gps, 2, Physics>: TemperatureGradient, Physics {
	TemperatureGradientKernel(const TemperatureGradient &base): TemperatureGradient(base) {}

	void move(int n)
	{
		gradient += 2 * n;
	}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = gradient;
		SIMD g0 = zeros(), g1 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t n = 0; n < nodes; ++n) {
				g0 = g0 + element.dND[gp][n][0] * element.temp[n];
				g1 = g1 + element.dND[gp][n][1] * element.temp[n];
			}
		}
		SIMD scale = load1(1. / gps);
		g0 = g0 * scale;
		g1 = g1 * scale;
		size_t size = std::min((size_t)SIMD::size, (size_t)(end - gradient));
		for (size_t s = 0; s < size; ++s) {
			out[2 * s + 0] = g0[s];
			out[2 * s + 1] = g1[s];
		}
		gradient += 2 * size;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct TemperatureGradientKernel<nodes, gps, 3, Physics>: TemperatureGradient, Physics {
	TemperatureGradientKernel(const TemperatureGradient &base): TemperatureGradient(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = gradient;
		SIMD g0 = zeros(), g1 = zeros(), g2 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t n = 0; n < nodes; ++n) {
				g0 = g0 + element.dND[gp][n][0] * element.temp[n];
				g1 = g1 + element.dND[gp][n][1] * element.temp[n];
				g2 = g2 + element.dND[gp][n][2] * element.temp[n];
			}
		}
		SIMD scale = load1(1. / gps);
		g0 = g0 * scale;
		g1 = g1 * scale;
		g2 = g2 * scale;
		size_t size = std::min((size_t)SIMD::size, (size_t)(end - gradient));
		for (size_t s = 0; s < size; ++s) {
			out[3 * s + 0] = g0[s];
			out[3 * s + 1] = g1[s];
			out[3 * s + 2] = g2[s];
		}
		gradient += 3 * size;
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_GRADIENT_H_ */
