
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_ADVECTION_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_ADVECTION_H_

#include "subkernels.h"

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t etype, class Physics> struct AdvectionSubkernel;

template <size_t nodes, size_t gps, size_t ndim, class Physics>
struct AdvectionSubkernel<nodes, gps, ndim, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Physics>: AdvectionKernel, Physics {
	AdvectionSubkernel(const AdvectionKernel &base): AdvectionKernel(base) {}

	void simd(typename Physics::Element &element)
	{

	}
};

template <size_t nodes, size_t gps, size_t ndim, class Physics>
struct AdvectionSubkernel<nodes, gps, ndim, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: AdvectionKernel, Physics {
	AdvectionSubkernel(const AdvectionKernel &base): AdvectionKernel(base) {}

	void simd(typename Physics::Element &element)
	{

	}
};

template <size_t nodes, size_t gps, size_t ndim, class Physics>
struct AdvectionSubkernel<nodes, gps, ndim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC, Physics>: AdvectionKernel, Physics {
	AdvectionSubkernel(const AdvectionKernel &base): AdvectionKernel(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD usq, besq, be[nodes];
			for (size_t d = 0; d < ndim; ++d) {
				SIMD u = element.ecf.advection[gp][d] * element.ecf.density[gp] * element.ecf.heatCapacity[gp];
				usq = usq + u * u;
				for (size_t n = 0; n < nodes; ++n) {
					be[n] = be[n] + u * element.dND[gp][n][d];
				}
			}
			for (size_t n = 0; n < nodes; ++n) {
				besq = besq + be[n] * be[n];
			}

			SIMD unorm = sqrt(usq), benorm = sqrt(besq);
			SIMD urnorm = positive_guarded_recip(unorm), bernorm = positive_guarded_recip(benorm);

			SIMD C05 = load1(.5);
			SIMD C1 = load1(1);
			SIMD C2 = load1(2);
			SIMD he = C2 * unorm * bernorm;
			SIMD rhe = benorm * C05 * urnorm;
			SIMD Pe = C2 * element.conductivity[gp] * rhe * urnorm;
			SIMD tau = max(SIMD(), C1 - Pe);
			SIMD adv = C05 * he * tau * urnorm;

			// element.conductivity[gp] = element.conductivity[gp] + element.ecf.sigma * he * unorm;

			for (size_t n = 0; n < nodes; ++n) {
				SIMD scale = element.det[gp] * load1(element.w[gp]) * (load1(element.N[gp][n]) + adv * be[n]);
				for (size_t m = 0; m < nodes; ++m) {
					SIMD res = load(out + (n * nodes + m) * SIMD::size);
					res = res + be[m] * scale;
					store(out + (n * nodes + m) * SIMD::size, res);
				}
			}
		}
		K += SIMD::size * nodes * nodes;
	}
};

template <size_t nodes, size_t gps, size_t ndim, class Physics>
struct AdvectionSubkernel<nodes, gps, ndim, HeatTransferElementType::ASYMMETRIC_GENERAL, Physics>: AdvectionKernel, Physics {
	AdvectionSubkernel(const AdvectionKernel &base): AdvectionKernel(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD usq, besq, be[nodes];
			for (size_t d = 0; d < ndim; ++d) {
				SIMD u = element.ecf.advection[gp][d] * element.ecf.density[gp] * element.ecf.heatCapacity[gp];
				usq = usq + u * u;
				for (size_t n = 0; n < nodes; ++n) {
					be[n] = be[n] + u * element.dND[gp][n][d];
				}
			}
			for (size_t n = 0; n < nodes; ++n) {
				besq = besq + be[n] * be[n];
			}

			SIMD unorm = sqrt(usq), benorm = sqrt(besq);
			SIMD urnorm = positive_guarded_recip(unorm), bernorm = positive_guarded_recip(benorm);

			SIMD C05 = load1(.5);
			SIMD C1 = load1(1);
			SIMD C2 = load1(2);
			SIMD he = C2 * unorm * bernorm;
			SIMD rhe = benorm * C05 * urnorm;
			SIMD Pe = C2 * element.conductivity[gp][0] * rhe * urnorm;
			SIMD tau = max(SIMD(), C1 - Pe);
			SIMD adv = C05 * he * tau * urnorm;

			for (size_t d = 0; d < ndim; ++d) {
				// element.conductivity[gp][d * ndim + d] = element.conductivity[gp][d * ndim + d] + element.ecf.sigma * he * unorm;
			}

			for (size_t n = 0; n < nodes; ++n) {
				SIMD scale = element.det[gp] * load1(element.w[gp]) * (load1(element.N[gp][n]) + adv * be[n]);
				for (size_t m = 0; m < nodes; ++m) {
					SIMD res = load(out + (n * nodes + m) * SIMD::size);
					res = res + be[m] * scale;
					store(out + (n * nodes + m) * SIMD::size, res);
				}
			}
		}
		K += SIMD::size * nodes * nodes;
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_ADVECTION_H_ */
