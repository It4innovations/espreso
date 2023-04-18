
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_

#include "subkernel.h"

namespace espreso {

struct Integration: SubKernel {
	const char* name() const { return "Integration"; }

	Integration()
	{
		isconst = false;
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE | Assembler::SOLUTION;
	}
};


template <size_t nodes, size_t gps, size_t ndim, size_t edim, class Physics> struct IntegrationKernel;

template <size_t nodes, size_t gps, class Physics>
struct IntegrationKernel<nodes, gps, 2, 2, Physics>: Integration, Physics {
	IntegrationKernel(const Integration &base): Integration(base) {}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD jacobian0 = zeros(), jacobian1 = zeros(), jacobian2 = zeros(), jacobian3 = zeros();

			for (size_t n = 0; n < nodes; ++n) {
				SIMD coordsX = element.coords[n][0];
				SIMD coordsY = element.coords[n][1];
				SIMD dNX = load1(element.dN[gp][n][0]);
				SIMD dNY = load1(element.dN[gp][n][1]);

				jacobian0 = jacobian0 + dNX * coordsX;
				jacobian1 = jacobian1 + dNX * coordsY;
				jacobian2 = jacobian2 + dNY * coordsX;
				jacobian3 = jacobian3 + dNY * coordsY;
			}

			element.det[gp] = jacobian0 * jacobian3 - jacobian1 * jacobian2;

			SIMD detJx = ones() / element.det[gp];
			SIMD inv0 =  detJx * jacobian3;
			SIMD inv1 = -detJx * jacobian1;
			SIMD inv2 = -detJx * jacobian2;
			SIMD inv3 =  detJx * jacobian0;

			for (size_t n = 0; n < nodes; ++n) {
				SIMD dNX = load1(element.dN[gp][n][0]);
				SIMD dNY = load1(element.dN[gp][n][1]);
				element.dND[gp][n][0] = inv0 * dNX + inv1 * dNY;
				element.dND[gp][n][1] = inv2 * dNX + inv3 * dNY;
			}
		}
	}
};

template <size_t nodes, size_t gps, class Physics>
struct IntegrationKernel<nodes, gps, 3, 3, Physics>: Integration, Physics {
	IntegrationKernel(const Integration &base): Integration(base) {}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD jacobian0 = zeros(), jacobian1 = zeros(), jacobian2 = zeros(), jacobian3 = zeros(), jacobian4 = zeros();
			SIMD jacobian5 = zeros(), jacobian6 = zeros(), jacobian7 = zeros(), jacobian8 = zeros();

			for (size_t n = 0; n < nodes; ++n) {
				SIMD coordsX = element.coords[n][0];
				SIMD coordsY = element.coords[n][1];
				SIMD coordsZ = element.coords[n][2];
				SIMD dNX = load1(element.dN[gp][n][0]);
				SIMD dNY = load1(element.dN[gp][n][1]);
				SIMD dNZ = load1(element.dN[gp][n][2]);

				jacobian0 = jacobian0 + dNX * coordsX;
				jacobian1 = jacobian1 + dNX * coordsY;
				jacobian2 = jacobian2 + dNX * coordsZ;
				jacobian3 = jacobian3 + dNY * coordsX;
				jacobian4 = jacobian4 + dNY * coordsY;
				jacobian5 = jacobian5 + dNY * coordsZ;
				jacobian6 = jacobian6 + dNZ * coordsX;
				jacobian7 = jacobian7 + dNZ * coordsY;
				jacobian8 = jacobian8 + dNZ * coordsZ;
			}

			element.det[gp] =
					+ jacobian0 * jacobian4 * jacobian8
					+ jacobian1 * jacobian5 * jacobian6
					+ jacobian2 * jacobian3 * jacobian7
					- jacobian2 * jacobian4 * jacobian6
					- jacobian1 * jacobian3 * jacobian8
					- jacobian0 * jacobian5 * jacobian7;

			SIMD detJx = ones() / element.det[gp];
			SIMD inv0 = detJx * ( jacobian8 * jacobian4 - jacobian7 * jacobian5);
			SIMD inv1 = detJx * (-jacobian8 * jacobian1 + jacobian7 * jacobian2);
			SIMD inv2 = detJx * ( jacobian5 * jacobian1 - jacobian4 * jacobian2);
			SIMD inv3 = detJx * (-jacobian8 * jacobian3 + jacobian6 * jacobian5);
			SIMD inv4 = detJx * ( jacobian8 * jacobian0 - jacobian6 * jacobian2);
			SIMD inv5 = detJx * (-jacobian5 * jacobian0 + jacobian3 * jacobian2);
			SIMD inv6 = detJx * ( jacobian7 * jacobian3 - jacobian6 * jacobian4);
			SIMD inv7 = detJx * (-jacobian7 * jacobian0 + jacobian6 * jacobian1);
			SIMD inv8 = detJx * ( jacobian4 * jacobian0 - jacobian3 * jacobian1);

			for (size_t n = 0; n < nodes; ++n) {
				SIMD dNX = load1(element.dN[gp][n][0]);
				SIMD dNY = load1(element.dN[gp][n][1]);
				SIMD dNZ = load1(element.dN[gp][n][2]);
				element.dND[gp][n][0] = inv0 * dNX + inv1 * dNY + inv2 * dNZ;
				element.dND[gp][n][1] = inv3 * dNX + inv4 * dNY + inv5 * dNZ;
				element.dND[gp][n][2] = inv6 * dNX + inv7 * dNY + inv8 * dNZ;
			}
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_ */
