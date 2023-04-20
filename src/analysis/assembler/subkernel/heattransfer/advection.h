
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_ADVECTION_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_ADVECTION_H_

#include "subkernels.h"

namespace espreso {

struct Advection: SubKernel {
	const char* name() const { return "Advection"; }

	ECFExpressionVector *expression;
	double *K, sigma;

	Advection()
	: expression(nullptr), K(nullptr), sigma(0)
	{
		isconst = false;
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE;
	}

	void activate(ECFExpressionVector *expression, double *K, double sigma)
	{
		this->expression = expression;
		this->K = K;
		this->sigma = sigma;
		if (this->expression) {
			isactive = 1;
		}
	}
};

//template <size_t nodes, size_t gps, size_t ndim, class Physics> struct AdvectionParams;
//
//template <size_t nodes, size_t gps, class Physics>
//struct AdvectionParams<nodes, gps, 2, Physics>: Physics {
//	void simd(typename Physics::Element &element)
//	{
//
//	}
//};
//
//template <size_t nodes, size_t gps, class Physics>
//struct AdvectionParams<nodes, gps, 3, Physics>: Physics {
//	void simd(typename Physics::Element &element)
//	{
//
//	}
//};
//
//template <size_t nodes, size_t gps, size_t ndim, enum ThermalConductivityConfiguration::MODEL model, class Physics> struct AdvectionKernel;
//
//template <size_t nodes, size_t gps, class Physics>
//struct AdvectionKernel<nodes, gps, 2, ThermalConductivityConfiguration::MODEL::DIAGONAL, Physics>: Advection, Physics {
//	AdvectionKernel(const Advection &base): Advection(base) {}
//
//	void simd(typename Physics::Element &element)
//	{
//		for (size_t gp = 0; gp < gps; ++gp) {
//			SIMD usq, besq, be[nodes];
//			for (size_t d = 0; d < ndim; ++d) {
//				SIMD u = element.ecf.advection[gp][d] * element.ecf.density[gp] * element.ecf.heatCapacity[gp];
//				usq = usq + u * u;
//				for (size_t n = 0; n < nodes; ++n) {
//					be[n] = be[n] + u * element.dND[gp][n][d];
//				}
//			}
//			for (size_t n = 0; n < nodes; ++n) {
//				besq = besq + be[n] * be[n];
//			}
//
//			SIMD unorm = sqrt(usq), benorm = sqrt(besq);
//			SIMD urnorm = positive_guarded_recip(unorm), bernorm = positive_guarded_recip(benorm);
//
//			SIMD C05 = load1(.5);
//			SIMD C1 = load1(1);
//			SIMD C2 = load1(2);
//			SIMD he = C2 * unorm * bernorm;
//			SIMD rhe = benorm * C05 * urnorm;
//			SIMD Pe = C2 * element.conductivity[gp] * rhe * urnorm;
//			SIMD tau = max(SIMD(), C1 - Pe);
//			SIMD adv = C05 * he * tau * urnorm;
//
//			// element.conductivity[gp] = element.conductivity[gp] + element.ecf.sigma * he * unorm;
//
//			for (size_t n = 0; n < nodes; ++n) {
//				SIMD scale = element.det[gp] * load1(element.w[gp]) * (load1(element.N[gp][n]) + adv * be[n]);
//				for (size_t m = 0; m < nodes; ++m) {
//					SIMD res = load(out + (n * nodes + m) * SIMD::size);
//					res = res + be[m] * scale;
//					store(out + (n * nodes + m) * SIMD::size, res);
//				}
//			}
//		}
//	}
//};
//
//template <size_t nodes, size_t gps, class Physics>
//struct AdvectionKernel<nodes, gps, 3, ThermalConductivityConfiguration::MODEL::DIAGONAL, Physics>: Advection, Physics {
//	AdvectionKernel(const Advection &base): Advection(base) {}
//
//	void simd(typename Physics::Element &element)
//	{
//
//	}
//};
//
//template <size_t nodes, size_t gps, size_t ndim, class Physics>
//struct AdvectionKernel<nodes, gps, ndim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC, Physics>: Advection, Physics {
//	AdvectionKernel(const Advection &base): Advection(base) {}
//
//	void simd(typename Physics::Element &element)
//	{
//		double * __restrict__ out = K;
//		for (size_t gp = 0; gp < gps; ++gp) {
//			SIMD usq, besq, be[nodes];
//			for (size_t d = 0; d < ndim; ++d) {
//				SIMD u = element.ecf.advection[gp][d] * element.ecf.density[gp] * element.ecf.heatCapacity[gp];
//				usq = usq + u * u;
//				for (size_t n = 0; n < nodes; ++n) {
//					be[n] = be[n] + u * element.dND[gp][n][d];
//				}
//			}
//			for (size_t n = 0; n < nodes; ++n) {
//				besq = besq + be[n] * be[n];
//			}
//
//			SIMD unorm = sqrt(usq), benorm = sqrt(besq);
//			SIMD urnorm = positive_guarded_recip(unorm), bernorm = positive_guarded_recip(benorm);
//
//			SIMD C05 = load1(.5);
//			SIMD C1 = load1(1);
//			SIMD C2 = load1(2);
//			SIMD he = C2 * unorm * bernorm;
//			SIMD rhe = benorm * C05 * urnorm;
//			SIMD Pe = C2 * element.conductivity[gp] * rhe * urnorm;
//			SIMD tau = max(SIMD(), C1 - Pe);
//			SIMD adv = C05 * he * tau * urnorm;
//
//			// element.conductivity[gp] = element.conductivity[gp] + element.ecf.sigma * he * unorm;
//
//			for (size_t n = 0; n < nodes; ++n) {
//				SIMD scale = element.det[gp] * load1(element.w[gp]) * (load1(element.N[gp][n]) + adv * be[n]);
//				for (size_t m = 0; m < nodes; ++m) {
//					SIMD res = load(out + (n * nodes + m) * SIMD::size);
//					res = res + be[m] * scale;
//					store(out + (n * nodes + m) * SIMD::size, res);
//				}
//			}
//		}
//		K += SIMD::size * nodes * nodes;
//	}
//};
//
//template <size_t nodes, size_t gps, size_t ndim, class Physics>
//struct AdvectionKernel<nodes, gps, ndim, HeatTransferElementType::ASYMMETRIC_GENERAL, Physics>: Advection, Physics {
//	AdvectionKernel(const Advection &base): Advection(base) {}
//
//	void simd(typename Physics::Element &element)
//	{
//		double * __restrict__ out = K;
//		for (size_t gp = 0; gp < gps; ++gp) {
//			SIMD usq, besq, be[nodes];
//			for (size_t d = 0; d < ndim; ++d) {
//				SIMD u = element.ecf.advection[gp][d] * element.ecf.density[gp] * element.ecf.heatCapacity[gp];
//				usq = usq + u * u;
//				for (size_t n = 0; n < nodes; ++n) {
//					be[n] = be[n] + u * element.dND[gp][n][d];
//				}
//			}
//			for (size_t n = 0; n < nodes; ++n) {
//				besq = besq + be[n] * be[n];
//			}
//
//			SIMD unorm = sqrt(usq), benorm = sqrt(besq);
//			SIMD urnorm = positive_guarded_recip(unorm), bernorm = positive_guarded_recip(benorm);
//
//			SIMD C05 = load1(.5);
//			SIMD C1 = load1(1);
//			SIMD C2 = load1(2);
//			SIMD he = C2 * unorm * bernorm;
//			SIMD rhe = benorm * C05 * urnorm;
//			SIMD Pe = C2 * element.conductivity[gp][0] * rhe * urnorm;
//			SIMD tau = max(SIMD(), C1 - Pe);
//			SIMD adv = C05 * he * tau * urnorm;
//
//			for (size_t d = 0; d < ndim; ++d) {
//				// element.conductivity[gp][d * ndim + d] = element.conductivity[gp][d * ndim + d] + element.ecf.sigma * he * unorm;
//			}
//
//			for (size_t n = 0; n < nodes; ++n) {
//				SIMD scale = element.det[gp] * load1(element.w[gp]) * (load1(element.N[gp][n]) + adv * be[n]);
//				for (size_t m = 0; m < nodes; ++m) {
//					SIMD res = load(out + (n * nodes + m) * SIMD::size);
//					res = res + be[m] * scale;
//					store(out + (n * nodes + m) * SIMD::size, res);
//				}
//			}
//		}
//		K += SIMD::size * nodes * nodes;
//	}
//};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_ADVECTION_H_ */
