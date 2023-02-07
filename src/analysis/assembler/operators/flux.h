
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_FLUX_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_FLUX_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/module/heattransfer.element.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

struct TemperatureFluxBase: ActionOperator {
	TemperatureFluxBase(size_t interval, NamedData *flux)
	: flux(flux->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin)
	{
		isconst = false;
		action = Action::SOLUTION;
	}

	double* flux;

	void move(int n)
	{
		flux += n;
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct TemperatureFlux;

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct TemperatureFlux<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Physics>: TemperatureFluxBase, Physics {
	using TemperatureFluxBase::TemperatureFluxBase;

	constexpr static double scale = 1. / gps;

	void sisd(typename Physics::Element &element)
	{
		flux[0] = flux[1] = 0;
		for (size_t gp = 0; gp < gps; ++gp) {
			double f[2] = { 0, 0 };
			for (size_t n = 0; n < nodes; ++n) {
				f[0] += element.dND[gp][n][0] * element.temp[n];
				f[1] += element.dND[gp][n][1] * element.temp[n];
			}
			flux[0] += f[0] * element.conductivity[gp];
			flux[1] += f[1] * element.conductivity[gp];
		}
		flux[0] *= scale;
		flux[1] *= scale;
		move(2);
	}

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
			f0 = f0 + fgp0 * element.conductivity[gp];
			f1 = f1 + fgp1 * element.conductivity[gp];
		}
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		for (size_t s = 0; s < SIMD::size; ++s) {
			out[2 * s + 0] = f0[s];
			out[2 * s + 1] = f1[s];
		}
		move(2 * SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		double * __restrict__ out = flux;
		SIMD f0 = zeros(), f1 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD fgp0 = zeros(), fgp1 = zeros();
			for (size_t n = 0; n < nodes; ++n) {
				fgp0 = fgp0 + element.dND[gp][n][0] * element.temp[n];
				fgp1 = fgp1 + element.dND[gp][n][1] * element.temp[n];
			}
			f0 = f0 + fgp0 * element.conductivity[gp];
			f1 = f1 + fgp1 * element.conductivity[gp];
		}
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		for (size_t s = 0; s <size; ++s) {
			out[2 * s + 0] = f0[s];
			out[2 * s + 1] = f1[s];
		}
		move(2 * size);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct TemperatureFlux<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Physics>: TemperatureFluxBase, Physics {
	using TemperatureFluxBase::TemperatureFluxBase;

	constexpr static double scale = 1. / gps;

	void sisd(typename Physics::Element &element)
	{
		flux[0] = flux[1] = flux[2] = 0;
		for (size_t gp = 0; gp < gps; ++gp) {
			double f[2] = { 0, 0 };
			for (size_t n = 0; n < nodes; ++n) {
				f[0] += element.dND[gp][n][0] * element.temp[n];
				f[1] += element.dND[gp][n][1] * element.temp[n];
				f[2] += element.dND[gp][n][2] * element.temp[n];
			}
			flux[0] += f[0] * element.conductivity[gp];
			flux[1] += f[1] * element.conductivity[gp];
			flux[2] += f[2] * element.conductivity[gp];
		}
		flux[0] *= scale;
		flux[1] *= scale;
		flux[2] *= scale;
		move(3);
	}

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
			f0 = f0 + fgp0 * element.conductivity[gp];
			f1 = f1 + fgp1 * element.conductivity[gp];
			f2 = f2 + fgp2 * element.conductivity[gp];
		}
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		f2 = f2 * sscale;
		for (size_t s = 0; s < SIMD::size; ++s) {
			out[3 * s + 0] = f0[s];
			out[3 * s + 1] = f1[s];
			out[3 * s + 2] = f2[s];
		}
		move(3 * SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
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
			f0 = f0 + fgp0 * element.conductivity[gp];
			f1 = f1 + fgp1 * element.conductivity[gp];
			f2 = f2 + fgp2 * element.conductivity[gp];
		}
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		f2 = f2 * sscale;
		for (size_t s = 0; s < size; ++s) {
			out[3 * s + 0] = f0[s];
			out[3 * s + 1] = f1[s];
			out[3 * s + 2] = f2[s];
		}
		move(3 * size);
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct TemperatureFlux<nodes, gps, 2, edim, etype, Physics>: TemperatureFluxBase, Physics {
	using TemperatureFluxBase::TemperatureFluxBase;

	constexpr static double scale = 1. / gps;

	void sisd(typename Physics::Element &element)
	{
		flux[0] = flux[1] = 0;
		for (size_t gp = 0; gp < gps; ++gp) {
			double f[2] = { 0, 0 };
			for (size_t n = 0; n < nodes; ++n) {
				f[0] += element.dND[gp][n][0] * element.temp[n];
				f[1] += element.dND[gp][n][1] * element.temp[n];
			}
			flux[0] += f[0] * element.conductivity[gp][0] + f[1] * element.conductivity[gp][1];
			flux[1] += f[0] * element.conductivity[gp][2] + f[1] * element.conductivity[gp][3];
		}
		flux[0] *= scale;
		flux[1] *= scale;
		move(2);
	}

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
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		for (size_t s = 0; s < SIMD::size; ++s) {
			out[2 * s + 0] = f0[s];
			out[2 * s + 1] = f1[s];
		}
		move(2 * SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
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
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		for (size_t s = 0; s <size; ++s) {
			out[2 * s + 0] = f0[s];
			out[2 * s + 1] = f1[s];
		}
		move(2 * size);
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct TemperatureFlux<nodes, gps, 3, edim, etype, Physics>: TemperatureFluxBase, Physics {
	using TemperatureFluxBase::TemperatureFluxBase;

	constexpr static double scale = 1. / gps;

	void sisd(typename Physics::Element &element)
	{
		flux[0] = flux[1] = flux[2] = 0;
		for (size_t gp = 0; gp < gps; ++gp) {
			double f[2] = { 0, 0 };
			for (size_t n = 0; n < nodes; ++n) {
				f[0] += element.dND[gp][n][0] * element.temp[n];
				f[1] += element.dND[gp][n][1] * element.temp[n];
				f[2] += element.dND[gp][n][2] * element.temp[n];
			}
			flux[0] += f[0] * element.conductivity[gp][0] + f[1] * element.conductivity[gp][1] + f[2] * element.conductivity[gp][2];
			flux[1] += f[0] * element.conductivity[gp][3] + f[1] * element.conductivity[gp][4] + f[2] * element.conductivity[gp][5];
			flux[2] += f[0] * element.conductivity[gp][6] + f[1] * element.conductivity[gp][7] + f[2] * element.conductivity[gp][8];
		}
		flux[0] *= scale;
		flux[1] *= scale;
		flux[2] *= scale;
		move(3);
	}

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
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		f2 = f2 * sscale;
		for (size_t s = 0; s < SIMD::size; ++s) {
			out[3 * s + 0] = f0[s];
			out[3 * s + 1] = f1[s];
			out[3 * s + 2] = f2[s];
		}
		move(3 * SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
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
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		f2 = f2 * sscale;
		for (size_t s = 0; s < size; ++s) {
			out[3 * s + 0] = f0[s];
			out[3 * s + 1] = f1[s];
			out[3 * s + 2] = f2[s];
		}
		move(3 * size);
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_FLUX_H_ */
