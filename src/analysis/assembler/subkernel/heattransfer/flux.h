
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_FLUX_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_FLUX_H_

#include "subkernel.h"

namespace espreso {

struct TemperatureFlux: SubKernel {
	const char* name() const { return "TemperatureFlux"; }

	double *flux, *end;

	TemperatureFlux()
	: flux(nullptr), end(nullptr)
	{
		isconst = false;
		action = Assembler::SOLUTION;
	}

	void activate(size_t interval, NamedData *flux)
	{
		this->flux = flux->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
		this->end = flux->data.data() + flux->data.size();
		isactive = 1;
	}
};

template <size_t nodes, size_t gps, size_t ndim, enum ThermalConductivityConfiguration::MODEL model, class Physics> struct TemperatureFluxKernel;

template <size_t nodes, size_t gps, class Physics>
struct TemperatureFluxKernel<nodes, gps, 2, ThermalConductivityConfiguration::MODEL::ISOTROPIC, Physics>: TemperatureFlux, Physics {
	TemperatureFluxKernel(const TemperatureFlux &base): TemperatureFlux(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp][0];
			f1 = f1 + fgp1 * element.conductivity[gp][0];
		}
		SIMD scale = load1(1. / gps);
		f0 = f0 * scale;
		f1 = f1 * scale;
		size_t size = std::min((size_t)SIMD::size, (size_t)(end - flux));
		for (size_t s = 0; s < size; ++s) {
			out[2 * s + 0] = f0[s];
			out[2 * s + 1] = f1[s];
		}
		flux += 2 * size;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct TemperatureFluxKernel<nodes, gps, 3, ThermalConductivityConfiguration::MODEL::ISOTROPIC, Physics>: TemperatureFlux, Physics {
	TemperatureFluxKernel(const TemperatureFlux &base): TemperatureFlux(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros(), f2 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros(), fgp2 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
				fgp2 = fgp2 + element.dND[gp][n][2] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp][0];
			f1 = f1 + fgp1 * element.conductivity[gp][0];
			f2 = f2 + fgp2 * element.conductivity[gp][0];
		}
		SIMD scale = load1(1. / gps);
		f0 = f0 * scale;
		f1 = f1 * scale;
		f2 = f2 * scale;
		size_t size = std::min((size_t)SIMD::size, (size_t)(end - flux));
		for (size_t s = 0; s < size; ++s) {
			out[3 * s + 0] = f0[s];
			out[3 * s + 1] = f1[s];
			out[3 * s + 2] = f2[s];
		}
		flux += 3 * size;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct TemperatureFluxKernel<nodes, gps, 2, ThermalConductivityConfiguration::MODEL::DIAGONAL, Physics>: TemperatureFlux, Physics {
	TemperatureFluxKernel(const TemperatureFlux &base): TemperatureFlux(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp][0];
			f1 = f1 + fgp1 * element.conductivity[gp][1];
		}
		SIMD scale = load1(1. / gps);
		f0 = f0 * scale;
		f1 = f1 * scale;
		size_t size = std::min((size_t)SIMD::size, (size_t)(end - flux));
		for (size_t s = 0; s < size; ++s) {
			out[2 * s + 0] = f0[s];
			out[2 * s + 1] = f1[s];
		}
		flux += 2 * size;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct TemperatureFluxKernel<nodes, gps, 3, ThermalConductivityConfiguration::MODEL::DIAGONAL, Physics>: TemperatureFlux, Physics {
	TemperatureFluxKernel(const TemperatureFlux &base): TemperatureFlux(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros(), f2 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros(), fgp2 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
				fgp2 = fgp2 + element.dND[gp][n][2] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp][0];
			f1 = f1 + fgp1 * element.conductivity[gp][1];
			f2 = f2 + fgp2 * element.conductivity[gp][2];
		}
		SIMD scale = load1(1. / gps);
		f0 = f0 * scale;
		f1 = f1 * scale;
		f2 = f2 * scale;
		size_t size = std::min((size_t)SIMD::size, (size_t)(end - flux));
		for (size_t s = 0; s < size; ++s) {
			out[3 * s + 0] = f0[s];
			out[3 * s + 1] = f1[s];
			out[3 * s + 2] = f2[s];
		}
		flux += 3 * size;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct TemperatureFluxKernel<nodes, gps, 2, ThermalConductivityConfiguration::MODEL::SYMMETRIC, Physics>: TemperatureFlux, Physics {
	TemperatureFluxKernel(const TemperatureFlux &base): TemperatureFlux(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp][0] + fgp1 * element.conductivity[gp][1];
			f1 = f1 + fgp0 * element.conductivity[gp][1] + fgp1 * element.conductivity[gp][2];
		}
		SIMD scale = load1(1. / gps);
		f0 = f0 * scale;
		f1 = f1 * scale;
		size_t size = std::min((size_t)SIMD::size, (size_t)(end - flux));
		for (size_t s = 0; s <size; ++s) {
			out[2 * s + 0] = f0[s];
			out[2 * s + 1] = f1[s];
		}
		flux += 2 * size;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct TemperatureFluxKernel<nodes, gps, 3, ThermalConductivityConfiguration::MODEL::SYMMETRIC, Physics>: TemperatureFlux, Physics {
	TemperatureFluxKernel(const TemperatureFlux &base): TemperatureFlux(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros(), f2 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros(), fgp2 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
				fgp2 = fgp2 + element.dND[gp][n][2] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp][0] + fgp1 * element.conductivity[gp][1] + fgp2 * element.conductivity[gp][2];
			f1 = f1 + fgp0 * element.conductivity[gp][1] + fgp1 * element.conductivity[gp][3] + fgp2 * element.conductivity[gp][4];
			f2 = f2 + fgp0 * element.conductivity[gp][2] + fgp1 * element.conductivity[gp][4] + fgp2 * element.conductivity[gp][5];
		}
		SIMD scale = load1(1. / gps);
		f0 = f0 * scale;
		f1 = f1 * scale;
		f2 = f2 * scale;
		size_t size = std::min((size_t)SIMD::size, (size_t)(end - flux));
		for (size_t s = 0; s < size; ++s) {
			out[3 * s + 0] = f0[s];
			out[3 * s + 1] = f1[s];
			out[3 * s + 2] = f2[s];
		}
		flux += 3 * size;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct TemperatureFluxKernel<nodes, gps, 2, ThermalConductivityConfiguration::MODEL::ANISOTROPIC, Physics>: TemperatureFlux, Physics {
	TemperatureFluxKernel(const TemperatureFlux &base): TemperatureFlux(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp][0] + fgp1 * element.conductivity[gp][1];
			f1 = f1 + fgp0 * element.conductivity[gp][2] + fgp1 * element.conductivity[gp][3];
		}
		SIMD scale = load1(1. / gps);
		f0 = f0 * scale;
		f1 = f1 * scale;
		size_t size = std::min((size_t)SIMD::size, (size_t)(end - flux));
		for (size_t s = 0; s <size; ++s) {
			out[2 * s + 0] = f0[s];
			out[2 * s + 1] = f1[s];
		}
		flux += 2 * size;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct TemperatureFluxKernel<nodes, gps, 3, ThermalConductivityConfiguration::MODEL::ANISOTROPIC, Physics>: TemperatureFlux, Physics {
	TemperatureFluxKernel(const TemperatureFlux &base): TemperatureFlux(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros(), f2 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros(), fgp2 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
				fgp2 = fgp2 + element.dND[gp][n][2] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp][0] + fgp1 * element.conductivity[gp][1] + fgp2 * element.conductivity[gp][2];
			f1 = f1 + fgp0 * element.conductivity[gp][3] + fgp1 * element.conductivity[gp][4] + fgp2 * element.conductivity[gp][5];
			f2 = f2 + fgp0 * element.conductivity[gp][6] + fgp1 * element.conductivity[gp][7] + fgp2 * element.conductivity[gp][8];
		}
		SIMD scale = load1(1. / gps);
		f0 = f0 * scale;
		f1 = f1 * scale;
		f2 = f2 * scale;
		size_t size = std::min((size_t)SIMD::size, (size_t)(end - flux));
		for (size_t s = 0; s < size; ++s) {
			out[3 * s + 0] = f0[s];
			out[3 * s + 1] = f1[s];
			out[3 * s + 2] = f2[s];
		}
		flux += 3 * size;
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_FLUX_H_ */
