
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_H_

#include "analysis/assembler/operator.h"
#include "math/simd/simd.h"

#include <cstdio>

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ElasticityIsotropicPlaneStrain: ActionOperator, Physics {

	ElasticityIsotropicPlaneStrain(size_t interval)
	{
		action = Action::ASSEMBLE;
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD ex = element.ecf.youngModulus[gp];
			SIMD mi = element.ecf.poissonRatio[gp];
			SIMD C1 = load1(1.);
			SIMD C2 = load1(2.);
			SIMD k = ex * (C1 - mi) / ((C1 + mi) * (C1 - C2 * mi));
			element.ecf.elasticity[gp][0] = k;
			element.ecf.elasticity[gp][1] = k * (mi / (C1 - mi));
			element.ecf.elasticity[gp][2] = zeros();
			element.ecf.elasticity[gp][3] = k * (mi / (C1 - mi));
			element.ecf.elasticity[gp][4] = k;
			element.ecf.elasticity[gp][5] = zeros();
			element.ecf.elasticity[gp][6] = zeros();
			element.ecf.elasticity[gp][7] = zeros();
			element.ecf.elasticity[gp][8] = k * ((C1 - C2 * mi) / (C2 * (C1 - mi)));
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ElasticityIsotropicPlaneStress: ActionOperator, Physics {

	ElasticityIsotropicPlaneStress(size_t interval)
	{
		action = Action::ASSEMBLE;
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD ex = element.ecf.youngModulus[gp];
			SIMD mi = element.ecf.poissonRatio[gp];
			SIMD C1 = load1(1.);
			SIMD C2 = load1(.5);
			SIMD k = ex / (C1 - mi * mi);
			element.ecf.elasticity[gp][0] = k;
			element.ecf.elasticity[gp][1] = k * mi;
			element.ecf.elasticity[gp][2] = zeros();
			element.ecf.elasticity[gp][3] = k * mi;
			element.ecf.elasticity[gp][4] = k;
			element.ecf.elasticity[gp][5] = zeros();
			element.ecf.elasticity[gp][6] = zeros();
			element.ecf.elasticity[gp][7] = zeros();
			element.ecf.elasticity[gp][8] = k * ((C1 -  mi) * C2);
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ElasticityIsotropicPlaneAxisymmetric: ActionOperator, Physics {

	ElasticityIsotropicPlaneAxisymmetric(size_t interval)
	{
		action = Action::ASSEMBLE;
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD ex = element.ecf.youngModulus[gp];
			SIMD mi = element.ecf.poissonRatio[gp];
			SIMD C05 = load1(.5);
			SIMD C1 = load1(1.);
			SIMD C2 = load1(2.);
			SIMD k = ex / ((C1 + mi) * (C1 - C2 * mi));
			element.ecf.elasticity[gp][ 0] = k * (C1 - mi);
			element.ecf.elasticity[gp][ 1] = k * mi;
			element.ecf.elasticity[gp][ 2] = k * mi;
			element.ecf.elasticity[gp][ 3] = zeros();

			element.ecf.elasticity[gp][ 4] = k * mi;
			element.ecf.elasticity[gp][ 5] = k * (C1 - mi);
			element.ecf.elasticity[gp][ 6] = k * mi;
			element.ecf.elasticity[gp][ 7] = zeros();

			element.ecf.elasticity[gp][ 8] = k * mi;
			element.ecf.elasticity[gp][ 9] = k * mi;
			element.ecf.elasticity[gp][10] = k * (C1 - mi);
			element.ecf.elasticity[gp][11] = zeros();

			element.ecf.elasticity[gp][12] = zeros();
			element.ecf.elasticity[gp][13] = zeros();
			element.ecf.elasticity[gp][14] = zeros();
			element.ecf.elasticity[gp][15] = k * (C1 - C2 * mi) * C05;
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ElasticityIsotropicVolume: ActionOperator, Physics {

	ElasticityIsotropicVolume(size_t interval)
	{
		action = Action::ASSEMBLE;
	}

	void simd(typename Physics::Element &element)
	{
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
			element.ecf.elasticity[gp][ 6] = ee * mi;
			element.ecf.elasticity[gp][ 7] = ee * (C1 - mi);
			element.ecf.elasticity[gp][ 8] = ee * mi;
			element.ecf.elasticity[gp][ 9] = zeros();
			element.ecf.elasticity[gp][10] = zeros();
			element.ecf.elasticity[gp][11] = zeros();
			element.ecf.elasticity[gp][12] = ee * mi;
			element.ecf.elasticity[gp][13] = ee * mi;
			element.ecf.elasticity[gp][14] = ee * (C1 - mi);
			element.ecf.elasticity[gp][15] = zeros();
			element.ecf.elasticity[gp][16] = zeros();
			element.ecf.elasticity[gp][17] = zeros();
			element.ecf.elasticity[gp][18] = zeros();
			element.ecf.elasticity[gp][19] = zeros();
			element.ecf.elasticity[gp][20] = zeros();
			element.ecf.elasticity[gp][21] = ee * (C05 - mi);
			element.ecf.elasticity[gp][22] = zeros();
			element.ecf.elasticity[gp][23] = zeros();
			element.ecf.elasticity[gp][24] = zeros();
			element.ecf.elasticity[gp][25] = zeros();
			element.ecf.elasticity[gp][26] = zeros();
			element.ecf.elasticity[gp][27] = zeros();
			element.ecf.elasticity[gp][28] = ee * (C05 - mi);
			element.ecf.elasticity[gp][29] = zeros();
			element.ecf.elasticity[gp][30] = zeros();
			element.ecf.elasticity[gp][31] = zeros();
			element.ecf.elasticity[gp][32] = zeros();
			element.ecf.elasticity[gp][33] = zeros();
			element.ecf.elasticity[gp][34] = zeros();
			element.ecf.elasticity[gp][35] = ee * (C05 - mi);
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ElasticityOrthotropicVolume: ActionOperator, Physics {

	ElasticityOrthotropicVolume(size_t interval)
	{
		action = Action::ASSEMBLE;
	}

	void simd(typename Physics::Element &element)
	{
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

			element.ecf.elasticity[gp][ 6] = dxy;
			element.ecf.elasticity[gp][ 7] = dyy;
			element.ecf.elasticity[gp][ 8] = dyz;
			element.ecf.elasticity[gp][ 9] = zeros();
			element.ecf.elasticity[gp][10] = zeros();
			element.ecf.elasticity[gp][11] = zeros();

			element.ecf.elasticity[gp][12] = dxz;
			element.ecf.elasticity[gp][13] = dyz;
			element.ecf.elasticity[gp][14] = dzz;
			element.ecf.elasticity[gp][15] = zeros();
			element.ecf.elasticity[gp][16] = zeros();
			element.ecf.elasticity[gp][17] = zeros();

			element.ecf.elasticity[gp][18] = zeros();
			element.ecf.elasticity[gp][19] = zeros();
			element.ecf.elasticity[gp][20] = zeros();
			element.ecf.elasticity[gp][21] = gx;
			element.ecf.elasticity[gp][22] = zeros();
			element.ecf.elasticity[gp][23] = zeros();

			element.ecf.elasticity[gp][24] = zeros();
			element.ecf.elasticity[gp][25] = zeros();
			element.ecf.elasticity[gp][26] = zeros();
			element.ecf.elasticity[gp][27] = zeros();
			element.ecf.elasticity[gp][28] = gy;
			element.ecf.elasticity[gp][29] = zeros();

			element.ecf.elasticity[gp][30] = zeros();
			element.ecf.elasticity[gp][31] = zeros();
			element.ecf.elasticity[gp][32] = zeros();
			element.ecf.elasticity[gp][33] = zeros();
			element.ecf.elasticity[gp][34] = zeros();
			element.ecf.elasticity[gp][35] = gz;
		}
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_H_ */
