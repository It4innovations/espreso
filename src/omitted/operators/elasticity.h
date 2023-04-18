
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_H_

#include <analysis/assembler/subkernel/operator.h>
#include "math/simd/simd.h"

#include <cstdio>

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ElasticityIsotropicPlaneStrain: ActionOperator, Physics {
	const char* name() const { return "ElasticityIsotropicPlaneStrain"; }

	ElasticityIsotropicPlaneStrain(size_t interval)
	{
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	void simd(typename Physics::Element &element)
	{
		// 0 1 2  0 1 2
		// 3 4 5    3 4
		// 6 7 8      5
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD ex = element.ecf.youngModulus[gp];
			SIMD mi = element.ecf.poissonRatio[gp];
			SIMD C1 = load1(1.);
			SIMD C2 = load1(2.);
			SIMD k = ex * (C1 - mi) / ((C1 + mi) * (C1 - C2 * mi));
			element.ecf.elasticity[gp][0] = k;
			element.ecf.elasticity[gp][1] = k * (mi / (C1 - mi));
			element.ecf.elasticity[gp][2] = zeros();
			element.ecf.elasticity[gp][3] = k;
			element.ecf.elasticity[gp][4] = zeros();
			element.ecf.elasticity[gp][5] = k * ((C1 - C2 * mi) / (C2 * (C1 - mi)));
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ElasticityIsotropicPlaneStress: ActionOperator, Physics {
	const char* name() const { return "ElasticityIsotropicPlaneStress"; }

	ElasticityIsotropicPlaneStress(size_t interval)
	{
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	void simd(typename Physics::Element &element)
	{
		// 0 1 2  0 1 2
		// 3 4 5    3 4
		// 6 7 8      5
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD ex = element.ecf.youngModulus[gp];
			SIMD mi = element.ecf.poissonRatio[gp];
			SIMD C1 = load1(1.);
			SIMD C2 = load1(.5);
			SIMD k = ex / (C1 - mi * mi);
			element.ecf.elasticity[gp][0] = k;
			element.ecf.elasticity[gp][1] = k * mi;
			element.ecf.elasticity[gp][2] = zeros();
			element.ecf.elasticity[gp][3] = k;
			element.ecf.elasticity[gp][4] = zeros();
			element.ecf.elasticity[gp][5] = k * ((C1 -  mi) * C2);
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ElasticityIsotropicPlaneAxisymmetric: ActionOperator, Physics {
	const char* name() const { return "ElasticityIsotropicPlaneAxisymmetric"; }

	ElasticityIsotropicPlaneAxisymmetric(size_t interval)
	{
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	void simd(typename Physics::Element &element)
	{
		//  0  1  2  3   0  1  2  3
		//  4  5  6  7      4  5  6
		//  8  9 10 11         7  8
		// 12 13 14 15            9
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD ex = element.ecf.youngModulus[gp];
			SIMD mi = element.ecf.poissonRatio[gp];
			SIMD C05 = load1(.5);
			SIMD C1 = load1(1.);
			SIMD C2 = load1(2.);
			SIMD k = ex / ((C1 + mi) * (C1 - C2 * mi));
			element.ecf.elasticity[gp][0] = k * (C1 - mi);
			element.ecf.elasticity[gp][1] = k * mi;
			element.ecf.elasticity[gp][2] = k * mi;
			element.ecf.elasticity[gp][3] = zeros();

			element.ecf.elasticity[gp][4] = k * (C1 - mi);
			element.ecf.elasticity[gp][5] = k * mi;
			element.ecf.elasticity[gp][6] = zeros();

			element.ecf.elasticity[gp][7] = k * (C1 - mi);
			element.ecf.elasticity[gp][8] = zeros();

			element.ecf.elasticity[gp][9] = k * (C1 - C2 * mi) * C05;
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ElasticityIsotropicVolume: ActionOperator, Physics {
	const char* name() const { return "ElasticityIsotropicVolume"; }

	ElasticityIsotropicVolume(size_t interval)
	{
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	void simd(typename Physics::Element &element)
	{
		//  0  1  2  3  4  5
		//  6  7  8  9 10 11
		// 12 13 14 15 16 17
		// 18 19 20 21 22 23
		// 24 25 26 27 28 29
		// 30 31 32 33 34 35
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD ex = element.ecf.youngModulus[gp][0];
			SIMD mi = element.ecf.poissonRatio[gp][0];
			SIMD C05 = load1(.5);
			SIMD C1 = load1(1.);
			SIMD C2 = load1(2.);
			SIMD ee = ex / ((C1 + mi) * (C1 - C2 * mi));
			element.ecf.elasticity[gp][ 0] = ee * (C1 - mi);
			element.ecf.elasticity[gp][ 1] = ee * mi;
			element.ecf.elasticity[gp][ 2] = ee * mi;
			element.ecf.elasticity[gp][ 3] = zeros();
			element.ecf.elasticity[gp][ 4] = zeros();
			element.ecf.elasticity[gp][ 5] = zeros();

			element.ecf.elasticity[gp][ 6] = ee * (C1 - mi);
			element.ecf.elasticity[gp][ 7] = ee * mi;
			element.ecf.elasticity[gp][ 8] = zeros();
			element.ecf.elasticity[gp][ 9] = zeros();
			element.ecf.elasticity[gp][10] = zeros();

			element.ecf.elasticity[gp][11] = ee * (C1 - mi);
			element.ecf.elasticity[gp][12] = zeros();
			element.ecf.elasticity[gp][13] = zeros();
			element.ecf.elasticity[gp][14] = zeros();

			element.ecf.elasticity[gp][15] = ee * (C05 - mi);
			element.ecf.elasticity[gp][16] = zeros();
			element.ecf.elasticity[gp][17] = zeros();

			element.ecf.elasticity[gp][18] = ee * (C05 - mi);
			element.ecf.elasticity[gp][19] = zeros();

			element.ecf.elasticity[gp][20] = ee * (C05 - mi);
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ElasticityOrthotropicVolume: ActionOperator, Physics {
	const char* name() const { return "ElasticityOrthotropicVolume"; }

	ElasticityOrthotropicVolume(size_t interval)
	{
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	void simd(typename Physics::Element &element)
	{
		//  0  1  2  3  4  5
		//  6  7  8  9 10 11
		// 12 13 14 15 16 17
		// 18 19 20 21 22 23
		// 24 25 26 27 28 29
		// 30 31 32 33 34 35
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD ex = element.ecf.youngModulus[gp][0];
			SIMD ey = element.ecf.youngModulus[gp][1];
			SIMD ez = element.ecf.youngModulus[gp][2];
			SIMD miXY = element.ecf.poissonRatio[gp][0];
			SIMD miXZ = element.ecf.poissonRatio[gp][1];
			SIMD miYZ = element.ecf.poissonRatio[gp][2];
			SIMD gx = element.ecf.shearModulus[gp][0];
			SIMD gy = element.ecf.shearModulus[gp][1];
			SIMD gz = element.ecf.shearModulus[gp][2];

			SIMD miYX = miXY * ey / ex;
			SIMD miZY = miYZ * ez / ey;
			SIMD miZX = miXZ * ex / ez;

			SIMD C1 = load1(1);
			SIMD ksi = C1 / (C1 - (miXY * miYX + miYZ * miZY + miXZ * miZX) - (miXY * miYZ * miZX + miYX * miZY * miXZ));

			SIMD dxx = ksi * ex * (C1 - miYZ * miZY);
			SIMD dxy = ksi * ey * (miXY + miXZ * miZY);
			SIMD dxz = ksi * ez * (miXZ + miYZ * miXY);
			SIMD dyy = ksi * ey * (C1 - miXZ * miZX);
			SIMD dyz = ksi * ez * (miYZ + miYX * miXZ);
			SIMD dzz = ksi * ez * (C1 - miYX * miXY);

			element.ecf.elasticity[gp][ 0] = dxx;
			element.ecf.elasticity[gp][ 1] = dxy;
			element.ecf.elasticity[gp][ 2] = dxz;
			element.ecf.elasticity[gp][ 3] = zeros();
			element.ecf.elasticity[gp][ 4] = zeros();
			element.ecf.elasticity[gp][ 5] = zeros();

			element.ecf.elasticity[gp][ 6] = dyy;
			element.ecf.elasticity[gp][ 7] = dyz;
			element.ecf.elasticity[gp][ 8] = zeros();
			element.ecf.elasticity[gp][ 9] = zeros();
			element.ecf.elasticity[gp][10] = zeros();

			element.ecf.elasticity[gp][11] = dzz;
			element.ecf.elasticity[gp][12] = zeros();
			element.ecf.elasticity[gp][13] = zeros();
			element.ecf.elasticity[gp][14] = zeros();

			element.ecf.elasticity[gp][15] = gx;
			element.ecf.elasticity[gp][16] = zeros();
			element.ecf.elasticity[gp][17] = zeros();

			element.ecf.elasticity[gp][18] = gy;
			element.ecf.elasticity[gp][19] = zeros();

			element.ecf.elasticity[gp][20] = gz;
		}
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_H_ */
