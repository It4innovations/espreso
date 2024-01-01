
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_STRESS_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_STRESS_H_

#include "analysis/assembler/general/element.h"
#include "analysis/assembler/general/subkernel.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nameddata.h"

namespace espreso {

struct Stress: SubKernel {
	const char* name() const { return "Stress"; }

	Stress()
	: principalStress(nullptr), componentStress(nullptr), vonMisesStress(nullptr), vonMisesStressEnd(nullptr)
	{
		isconst = false;
		action = SubKernel::SOLUTION;
	}

	double* principalStress, *componentStress, *vonMisesStress, *vonMisesStressEnd;

	void activate(size_t interval, NamedData *principalStress, NamedData *componentStress, NamedData *vonMisesStress)
	{
		this->principalStress = principalStress->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
		this->componentStress = componentStress->data.data() + 2 * info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
		this->vonMisesStress = vonMisesStress->data.data() + info::mesh->elements->eintervals[interval].begin;
		this->vonMisesStressEnd = vonMisesStress->data.data() + vonMisesStress->data.size();
		isactive = 1;
	}
};

template <size_t nodes, size_t gps, size_t ndim> struct StressKernel;

template <size_t nodes, size_t gps>
struct StressKernel<nodes, gps, 2>: Stress {
	StressKernel(const Stress &base): Stress(base) {}

	template <typename Element>
	void simd(Element &element)
	{
		// TODO
	}
};

template <size_t nodes, size_t gps>
struct StressKernel<nodes, gps, 3>: Stress {
	StressKernel(const Stress &base): Stress(base) {}

	const double rgps = 1.0 / gps;

	template <typename Element>
	void simd(Element &element)
	{
		size_t size = std::min((size_t)SIMD::size, (size_t)(vonMisesStressEnd - vonMisesStress));

		SIMD scale = load1(rgps);
		SIMD CuB0 = element.sigma[0] * scale;
		SIMD CuB1 = element.sigma[1] * scale;
		SIMD CuB2 = element.sigma[2] * scale;
		SIMD CuB3 = element.sigma[3] * scale;
		SIMD CuB4 = element.sigma[4] * scale;
		SIMD CuB5 = element.sigma[5] * scale;

		double * __restrict__ component = componentStress;
		for (size_t s = 0; s < size; ++s) {
			component[6 * s + 0] = CuB0[s];
			component[6 * s + 1] = CuB1[s];
			component[6 * s + 2] = CuB2[s];
			component[6 * s + 3] = CuB3[s];
			component[6 * s + 4] = CuB4[s];
			component[6 * s + 5] = CuB5[s];
		}

		// https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
		// p1 = A(1,2)^2 + A(1,3)^2 + A(2,3)^2
		// if (p1 == 0)
		//    % A is diagonal.
		//    eig1 = A(1,1)
		//    eig2 = A(2,2)
		//    eig3 = A(3,3)
		// else
		//    q = trace(A)/3               % trace(A) is the sum of all diagonal values
		//    p2 = (A(1,1) - q)^2 + (A(2,2) - q)^2 + (A(3,3) - q)^2 + 2 * p1
		//    p = sqrt(p2 / 6)
		//    B = (1 / p) * (A - q * I)    % I is the identity matrix
		//    r = det(B) / 2
		//
		//    % In exact arithmetic for a symmetric matrix  -1 <= r <= 1
		//    % but computation error can leave it slightly outside this range.
		//    if (r <= -1)
		//       phi = pi / 3
		//    elseif (r >= 1)
		//       phi = 0
		//    else
		//       phi = acos(r) / 3
		//    end
		//
		//    % the eigenvalues satisfy eig3 <= eig2 <= eig1
		//    eig1 = q + 2 * p * cos(phi)
		//    eig3 = q + 2 * p * cos(phi + (2*pi/3))
		//    eig2 = 3 * q - eig1 - eig3     % since trace(A) = eig1 + eig2 + eig3
		// end



		// CuB0 CuB3 CuB5
		// CuB3 CuB1 CuB4
		// CuB5 cUB4 CuB2
		SIMD p1 = CuB3 * CuB3 + CuB4 * CuB4 + CuB5 * CuB5;
		SIMD q  = load1(1. / 3) * (CuB0 + CuB1 + CuB2);
		SIMD p2 = load1(1. / 6) * ((CuB0 - q) * (CuB0 - q) + (CuB1 - q) * (CuB1 - q) + (CuB2 - q) * (CuB2 - q) + load1(2) * p1);
		SIMD p  = sqrt(p2);
		SIMD rp = rsqrt14(p2);
		SIMD B0 = rp * (CuB0 - q);
		SIMD B1 = rp * (CuB1 - q);
		SIMD B2 = rp * (CuB2 - q);
		SIMD B3 = rp * CuB3;
		SIMD B4 = rp * CuB4;
		SIMD B5 = rp * CuB5;

		SIMD r = load1(1. / 2) * (B0 * B1 * B2 + load1(2) * B3 * B4 * B5 - B5 * B5 * B1 - B3 * B3 * B2 - B4 * B4 * B0);
		SIMD phi = load1(1. / 3) * acos(r);

		SIMD e1 = q + load1(2) * p * cos(phi);
		SIMD e3 = q + load1(2) * p * cos(phi + load1(2) * load1(M_PI) * load1(1. / 3));
		SIMD e2 = load1(3) * q - e1 - e3;

		double * __restrict__ principal = principalStress;
		for (size_t s = 0; s < size; ++s) {
			principal[3 * s + 0] = e1[s];
			principal[3 * s + 1] = e2[s];
			principal[3 * s + 2] = e3[s];
		}

		SIMD vm = (e1 - e2) * (e1 - e2) + (e1 - e3) * (e1 - e3) + (e2 - e3) * (e2 - e3);
		vm = sqrt(load1(1. / 2) * vm);
		double * __restrict__ vonMises = vonMisesStress;
		for (size_t s = 0; s < size; ++s) {
			vonMises[s] = vm[s];
		}

		principalStress += 3 * SIMD::size;
		componentStress += 6 * SIMD::size;
		vonMisesStress  += 1 * SIMD::size;
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_STRESS_H_ */
