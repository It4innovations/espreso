
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_STRESS_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_STRESS_H_

#include "subkernel.h"
#include "mesh/store/elementstore.h"
#include "esinfo/meshinfo.h"

namespace espreso {

struct Stress: SubKernel {
	const char* name() const { return "StressBase"; }

	Stress()
	: principalStress(nullptr), componentStress(nullptr), vonMisesStress(nullptr), vonMisesStressEnd(nullptr)
	{
		isconst = false;
		action = Assembler::SOLUTION;
	}

	double* principalStress, *componentStress, *vonMisesStress, *vonMisesStressEnd;

	void activate(size_t interval, NamedData *principalStress, NamedData *componentStress, NamedData *vonMisesStress)
	{
		this->principalStress = principalStress->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
		this->componentStress = componentStress->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
		this->vonMisesStress = vonMisesStress->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
		this->vonMisesStressEnd = vonMisesStress->data.data() + vonMisesStress->data.size();
		isactive = 1;
	}
};

template <size_t nodes, size_t gps, size_t ndim, ElasticityModel model, class Physics> struct StressKernel;

template <size_t nodes, size_t gps, ElasticityModel model, class Physics>
struct StressKernel<nodes, gps, 2, model, Physics>: Stress, Physics {
	StressKernel(const Stress &base): Stress(base) {}

	void simd(typename Physics::Element &element)
	{

	}
};

template <size_t nodes, class Physics>
void uB(typename Physics::Element &element, const size_t &gp, SIMD &uB0, SIMD &uB1, SIMD &uB2, SIMD &uB3, SIMD &uB4, SIMD &uB5)
{
	for (size_t n = 0; n < nodes; ++n) {
		SIMD dx = element.displacement[n][0];
		SIMD dy = element.displacement[n][1];
		SIMD dz = element.displacement[n][2];
		SIMD dNDx = element.dND[gp][n][0];
		SIMD dNDy = element.dND[gp][n][1];
		SIMD dNDz = element.dND[gp][n][2];

		uB0 = uB0 + dNDx * dx;
		uB1 = uB1 + dNDy * dy;
		uB2 = uB2 + dNDz * dz;
		uB3 = uB3 + dNDy * dx + dNDx * dy;
		uB4 = uB4 + dNDz * dy + dNDy * dz;
		uB5 = uB5 + dNDz * dx + dNDx * dz;
	}
}

template <size_t nodes, size_t gps, size_t ndim, ElasticityModel model, class Physics> struct CuB;

template <size_t nodes, size_t gps, class Physics>
struct CuB<nodes, gps, 3, ElasticityModel::ISOTROPIC, Physics> {
	// C
	// 0 1 1 _ _ _
	//   0 1 _ _ _
	//     0 _ _ _
	//       2 _ _
	//         2 _
	//           2
	static void simd(typename Physics::Element &element, SIMD &CuB0, SIMD &CuB1, SIMD &CuB2, SIMD &CuB3, SIMD &CuB4, SIMD &CuB5)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD c0 = element.elasticity[gp][0];
			SIMD c1 = element.elasticity[gp][1];
			SIMD c2 = element.elasticity[gp][2];

			SIMD uB0, uB1, uB2, uB3, uB4, uB5;
			uB<nodes, Physics>(element, gp, uB0, uB1, uB2, uB3, uB4, uB5);

			CuB0 = CuB0 + uB0 * c0 + uB1 * c1 + uB2 * c1;
			CuB1 = CuB1 + uB0 * c1 + uB1 * c0 + uB2 * c1;
			CuB2 = CuB2 + uB0 * c1 + uB1 * c1 + uB2 * c0;
			CuB3 = CuB3 + uB3 * c2;
			CuB4 = CuB4 + uB4 * c2;
			CuB5 = CuB5 + uB5 * c2;
		}
	}
};

template <size_t nodes, size_t gps, class Physics>
struct CuB<nodes, gps, 3, ElasticityModel::ORTHOTROPIC, Physics> {
	// C
	// 0 1 2 _ _ _
	//   3 4 _ _ _
	//     5 _ _ _
	//       6 _ _
	//         7 _
	//           8
	static void simd(typename Physics::Element &element, SIMD &CuB0, SIMD &CuB1, SIMD &CuB2, SIMD &CuB3, SIMD &CuB4, SIMD &CuB5)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD c00 = element.elasticity[gp][0];
			SIMD c01 = element.elasticity[gp][1], c11 = element.elasticity[gp][3];
			SIMD c02 = element.elasticity[gp][2], c12 = element.elasticity[gp][4], c22 = element.elasticity[gp][5];
			SIMD c33 = element.elasticity[gp][6];
			SIMD c44 = element.elasticity[gp][7];
			SIMD c55 = element.elasticity[gp][8];

			SIMD uB0, uB1, uB2, uB3, uB4, uB5;
			uB<nodes, Physics>(element, gp, uB0, uB1, uB2, uB3, uB4, uB5);

			CuB0 = CuB0 + uB0 * c00 + uB1 * c01 + uB2 * c02;
			CuB1 = CuB1 + uB0 * c01 + uB1 * c11 + uB2 * c12;
			CuB2 = CuB2 + uB0 * c02 + uB1 * c12 + uB2 * c22;
			CuB3 = CuB3 + uB3 * c33;
			CuB4 = CuB4 + uB4 * c44;
			CuB5 = CuB5 + uB5 * c55;
		}
	}
};

template <size_t nodes, size_t gps, class Physics>
struct CuB<nodes, gps, 3, ElasticityModel::SYMMETRIC, Physics> {
	static void simd(typename Physics::Element &element, SIMD &CuB0, SIMD &CuB1, SIMD &CuB2, SIMD &CuB3, SIMD &CuB4, SIMD &CuB5)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD c00 = element.elasticity[gp][0];
			SIMD c01 = element.elasticity[gp][1], c11 = element.elasticity[gp][ 6];
			SIMD c02 = element.elasticity[gp][2], c12 = element.elasticity[gp][ 7], c22 = element.elasticity[gp][11];
			SIMD c03 = element.elasticity[gp][3], c13 = element.elasticity[gp][ 8], c23 = element.elasticity[gp][12], c33 = element.elasticity[gp][15];
			SIMD c04 = element.elasticity[gp][4], c14 = element.elasticity[gp][ 9], c24 = element.elasticity[gp][13], c34 = element.elasticity[gp][16], c44 = element.elasticity[gp][18];
			SIMD c05 = element.elasticity[gp][5], c15 = element.elasticity[gp][10], c25 = element.elasticity[gp][14], c35 = element.elasticity[gp][17], c45 = element.elasticity[gp][19], c55 = element.elasticity[gp][20];

			SIMD uB0, uB1, uB2, uB3, uB4, uB5;
			uB<nodes, Physics>(element, gp, uB0, uB1, uB2, uB3, uB4, uB5);

			CuB0 = CuB0 + uB0 * c00 + uB1 * c01 + uB2 * c02 + uB3 * c03 + uB4 * c04 + uB5 * c05;
			CuB1 = CuB1 + uB0 * c01 + uB1 * c11 + uB2 * c12 + uB3 * c13 + uB4 * c14 + uB5 * c15;
			CuB2 = CuB2 + uB0 * c02 + uB1 * c12 + uB2 * c22 + uB3 * c23 + uB4 * c24 + uB5 * c25;
			CuB3 = CuB3 + uB0 * c03 + uB1 * c13 + uB2 * c23 + uB3 * c33 + uB4 * c34 + uB5 * c35;
			CuB4 = CuB4 + uB0 * c04 + uB1 * c14 + uB2 * c24 + uB3 * c34 + uB4 * c44 + uB5 * c45;
			CuB5 = CuB5 + uB0 * c05 + uB1 * c15 + uB2 * c25 + uB3 * c35 + uB4 * c45 + uB5 * c55;
		}
	}
};

template <size_t nodes, size_t gps, class Physics>
struct CuB<nodes, gps, 3, ElasticityModel::ANISOTROPIC, Physics> {
	static void simd(typename Physics::Element &element, SIMD &CuB0, SIMD &CuB1, SIMD &CuB2, SIMD &CuB3, SIMD &CuB4, SIMD &CuB5)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD c00 = element.elasticity[gp][0], c10 = element.elasticity[gp][ 6], c20 = element.elasticity[gp][12], c30 = element.elasticity[gp][18], c40 = element.elasticity[gp][24], c50 = element.elasticity[gp][30];
			SIMD c01 = element.elasticity[gp][1], c11 = element.elasticity[gp][ 7], c21 = element.elasticity[gp][13], c31 = element.elasticity[gp][19], c41 = element.elasticity[gp][25], c51 = element.elasticity[gp][31];
			SIMD c02 = element.elasticity[gp][2], c12 = element.elasticity[gp][ 8], c22 = element.elasticity[gp][14], c32 = element.elasticity[gp][20], c42 = element.elasticity[gp][26], c52 = element.elasticity[gp][32];
			SIMD c03 = element.elasticity[gp][3], c13 = element.elasticity[gp][ 9], c23 = element.elasticity[gp][15], c33 = element.elasticity[gp][21], c43 = element.elasticity[gp][27], c53 = element.elasticity[gp][33];
			SIMD c04 = element.elasticity[gp][4], c14 = element.elasticity[gp][10], c24 = element.elasticity[gp][16], c34 = element.elasticity[gp][22], c44 = element.elasticity[gp][28], c54 = element.elasticity[gp][34];
			SIMD c05 = element.elasticity[gp][5], c15 = element.elasticity[gp][11], c25 = element.elasticity[gp][17], c35 = element.elasticity[gp][23], c45 = element.elasticity[gp][29], c55 = element.elasticity[gp][35];

			SIMD uB0, uB1, uB2, uB3, uB4, uB5;
			uB<nodes, Physics>(element, gp, uB0, uB1, uB2, uB3, uB4, uB5);

			CuB0 = CuB0 + uB0 * c00 + uB1 * c10 + uB2 * c20 + uB3 * c30 + uB4 * c40 + uB5 * c50;
			CuB1 = CuB1 + uB0 * c01 + uB1 * c11 + uB2 * c21 + uB3 * c31 + uB4 * c41 + uB5 * c51;
			CuB2 = CuB2 + uB0 * c02 + uB1 * c12 + uB2 * c22 + uB3 * c32 + uB4 * c42 + uB5 * c52;
			CuB3 = CuB3 + uB0 * c03 + uB1 * c13 + uB2 * c23 + uB3 * c33 + uB4 * c43 + uB5 * c53;
			CuB4 = CuB4 + uB0 * c04 + uB1 * c14 + uB2 * c24 + uB3 * c34 + uB4 * c44 + uB5 * c54;
			CuB5 = CuB5 + uB0 * c05 + uB1 * c15 + uB2 * c25 + uB3 * c35 + uB4 * c45 + uB5 * c55;
		}
	}
};

template <size_t nodes, size_t gps, ElasticityModel model, class Physics>
struct StressKernel<nodes, gps, 3, model, Physics>: Stress, Physics {
	StressKernel(const Stress &base): Stress(base) {}

	const double rgps = 1.0 / gps;

	// B = dX  0  0 dY  0 dZ
	//      0 dY  0 dX dZ  0
	//      0  0 dZ  0 dY dX


	void simd(typename Physics::Element &element)
	{
		size_t size = std::min((size_t)SIMD::size, (size_t)(vonMisesStressEnd - vonMisesStress));
		SIMD CuB0, CuB1, CuB2, CuB3, CuB4, CuB5;
		CuB<nodes, gps, 3, model, Physics>::simd(element, CuB0, CuB1, CuB2, CuB3, CuB4, CuB5);

		SIMD scale = load1(rgps);
		CuB0 = CuB0 * scale;
		CuB1 = CuB1 * scale;
		CuB2 = CuB2 * scale;
		CuB3 = CuB3 * scale;
		CuB4 = CuB4 * scale;
		CuB5 = CuB5 * scale;

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

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_STRESS_H_ */
