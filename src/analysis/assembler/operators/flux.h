
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_FLUX_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_FLUX_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/module/heattransfer.element.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

struct TemperatureFluxBase: ActionOperator {
	TemperatureFluxBase(int interval, NamedData *flux)
	: flux(flux->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin)
	{
		action = Action::SOLUTION;
	}

	double* flux;

	void move(int n)
	{
		flux += n * info::mesh->dimension;
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
			flux[0] += f[0] * element.conductivity[gp][0];
			flux[1] += f[1] * element.conductivity[gp][0];
		}
		flux[0] *= scale;
		flux[1] *= scale;
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
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
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		for (size_t s = 0; s < SIMD::size; ++s) {
			flux[2 * s + 0] = f0[s];
			flux[2 * s + 1] = f1[s];
		}
		move(SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
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
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		for (size_t s = 0; s <size; ++s) {
			flux[2 * s + 0] = f0[s];
			flux[2 * s + 1] = f1[s];
		}
		move(size);
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
			flux[0] += f[0] * element.conductivity[gp][0];
			flux[1] += f[1] * element.conductivity[gp][0];
			flux[2] += f[2] * element.conductivity[gp][0];
		}
		flux[0] *= scale;
		flux[1] *= scale;
		flux[2] *= scale;
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
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
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		f2 = f2 * sscale;
		for (size_t s = 0; s < SIMD::size; ++s) {
			flux[3 * s + 0] = f0[s];
			flux[3 * s + 1] = f1[s];
			flux[3 * s + 2] = f2[s];
		}
		move(SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
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
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		f2 = f2 * sscale;
		for (size_t s = 0; s < size; ++s) {
			flux[3 * s + 0] = f0[s];
			flux[3 * s + 1] = f1[s];
			flux[3 * s + 2] = f2[s];
		}
		move(size);
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
			flux[0] += f[0] * element.conductivity[gp][0];
			flux[1] += f[1] * element.conductivity[gp][0];
		}
		flux[0] *= scale;
		flux[1] *= scale;
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
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
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		for (size_t s = 0; s < SIMD::size; ++s) {
			flux[2 * s + 0] = f0[s];
			flux[2 * s + 1] = f1[s];
		}
		move(SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
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
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		for (size_t s = 0; s <size; ++s) {
			flux[2 * s + 0] = f0[s];
			flux[2 * s + 1] = f1[s];
		}
		move(size);
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
			flux[0] += f[0] * element.conductivity[gp][0];
			flux[1] += f[1] * element.conductivity[gp][0];
			flux[2] += f[2] * element.conductivity[gp][0];
		}
		flux[0] *= scale;
		flux[1] *= scale;
		flux[2] *= scale;
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
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
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		f2 = f2 * sscale;
		for (size_t s = 0; s < SIMD::size; ++s) {
			flux[3 * s + 0] = f0[s];
			flux[3 * s + 1] = f1[s];
			flux[3 * s + 2] = f2[s];
		}
		move(SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
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
		SIMD sscale; sscale.fill(scale);
		f0 = f0 * sscale;
		f1 = f1 * sscale;
		f2 = f2 * sscale;
		for (size_t s = 0; s < size; ++s) {
			flux[3 * s + 0] = f0[s];
			flux[3 * s + 1] = f1[s];
			flux[3 * s + 2] = f2[s];
		}
		move(size);
	}
};

}

//#include "analysis/assembler/operator.h"
//#include "analysis/assembler/parameter.h"
//#include "analysis/assembler/math.hpp"
//#include "esinfo/meshinfo.h"
//#include "mesh/store/elementstore.h"
//
//namespace espreso {
//
//struct OutputFlux: public ActionOperator {
//	OutputFlux(int interval, const ParameterData &dND, const ParameterData &temperature, const ParameterData &conductivity, NamedData *gradient)
//	: dND(dND, interval),
//	  temp(temperature, interval),
//	  conductivity(conductivity, interval),
//	  flux(gradient->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin)
//	{
//
//	}
//
//	InputParameterIterator dND, temp, conductivity;
//	double* flux;
//
//	void operator++()
//	{
//		++dND; ++temp; ++conductivity;
//		flux += info::mesh->dimension;
//	}
//
//	void move(int n)
//	{
//		dND += n; temp += n; conductivity += n;
//		flux += n * info::mesh->dimension;
//	}
//};
//
//template<size_t nodes, size_t gps>
//struct OutputFluxIsotropic2D: public OutputFlux {
//	using OutputFlux::OutputFlux;
//
//	void operator()()
//	{
//		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			ADDM2NMN1<nodes>(conductivity[gpindex] / gps, dND.data + 2 * nodes * gpindex, temp.data, flux);
//		}
//	}
//};
//
//template<size_t nodes, size_t gps>
//struct OutputFluxIsotropic3D: public OutputFlux {
//	using OutputFlux::OutputFlux;
//
//	void operator()()
//	{
//		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			ADDM2NMN1<nodes>(conductivity[gpindex] / gps, dND.data + 3 * nodes * gpindex, temp.data, flux);
//		}
//	}
//};
//
//template<size_t nodes, size_t gps>
//struct OutputFlux2D: public OutputFlux {
//	using OutputFlux::OutputFlux;
//
//	void operator()()
//	{
//		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			ADDM22M2NMN1<nodes>(1. / gps, conductivity.data + 4 * gpindex, dND.data + 2 * nodes * gpindex, temp.data, flux);
//		}
//	}
//};
//
//template<size_t nodes, size_t gps>
//struct OutputFlux3D: public OutputFlux {
//	using OutputFlux::OutputFlux;
//
//	void operator()()
//	{
//		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			ADDM33M3NMN1<nodes>(1. / gps, conductivity.data + 9 * gpindex, dND.data + 3 * nodes * gpindex, temp.data, flux);
//		}
//	}
//};
//
//}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_FLUX_H_ */
