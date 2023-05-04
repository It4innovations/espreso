
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_ADVECTION_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_ADVECTION_H_

#include "subkernel.h"

namespace espreso {

struct Advection: SubKernel {
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
			this->isactive = 1;
		}
	}
};


template <size_t ndim, enum ThermalConductivityConfiguration::MODEL model, class Physics> struct UpdateConductivity;

template <size_t ndim, class Physics> struct UpdateConductivity<ndim, ThermalConductivityConfiguration::MODEL::ISOTROPIC, Physics> {
	static void update(typename Physics::Element &element, size_t gp, SIMD stabilization)
	{
		element.conductivity[gp][0] = element.conductivity[gp][0] + stabilization;
	}
};

template <size_t ndim, class Physics> struct UpdateConductivity<ndim, ThermalConductivityConfiguration::MODEL::DIAGONAL, Physics> {
	static void update(typename Physics::Element &element, size_t gp, SIMD stabilization)
	{
		for (size_t d = 0; d < ndim; ++d) {
			element.conductivity[gp][d] = element.conductivity[gp][d] + stabilization;
		}
	}
};

template <size_t ndim, class Physics> struct UpdateConductivity<ndim, ThermalConductivityConfiguration::MODEL::SYMMETRIC, Physics> {
	static void update(typename Physics::Element &element, size_t gp, SIMD stabilization)
	{
		for (size_t d = 0, i = 0; d < ndim; i += (ndim - d++)) {
			element.conductivity[gp][i] = element.conductivity[gp][i] + stabilization;
		}
	}
};

template <size_t ndim, class Physics> struct UpdateConductivity<ndim, ThermalConductivityConfiguration::MODEL::ANISOTROPIC, Physics> {
	static void update(typename Physics::Element &element, size_t gp, SIMD stabilization)
	{
		for (size_t d = 0, i = 0; d < ndim; i += ++d + ndim) {
			element.conductivity[gp][i] = element.conductivity[gp][i] + stabilization;
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, enum ThermalConductivityConfiguration::MODEL model, class Physics> struct AdvectionKernel: Advection, Physics {
	AdvectionKernel(const Advection &base): Advection(base) {}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD usq, besq;
			for (size_t n = 0; n < nodes; ++n) {
				element.be[gp][n] = zeros();
			}
			for (size_t d = 0; d < ndim; ++d) {
				SIMD u = element.ecf.advection[gp][d] * element.ecf.density[gp] * element.ecf.heatCapacity[gp];
				usq = usq + u * u;
				for (size_t n = 0; n < nodes; ++n) {
					element.be[gp][n] = element.be[gp][n] + u * element.dND[gp][n][d];
				}
			}
			for (size_t n = 0; n < nodes; ++n) {
				besq = besq + element.be[gp][n] * element.be[gp][n];
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
			element.advection[gp] = C05 * he * tau * urnorm;

			UpdateConductivity<ndim, model, Physics>::update(element, gp, load1(sigma) * he * unorm);
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, enum ThermalConductivityConfiguration::MODEL model, class Physics> struct AdvectionMatrix: Advection, Physics {
	AdvectionMatrix(const Advection &base): Advection(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t n = 0; n < nodes; ++n) {
				SIMD scale = element.det[gp] * load1(element.w[gp]) * (load1(element.N[gp][n]) + element.advection[gp] * element.be[gp][n]);
				for (size_t m = 0; m < nodes; ++m) {
					SIMD res = load(out + (n * nodes + m) * SIMD::size);
					res = res + element.be[gp][m] * scale;
					store(out + (n * nodes + m) * SIMD::size, res);
				}
			}
		}
		K += SIMD::size * nodes * nodes;
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_ADVECTION_H_ */
