
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct Integration;

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct Integration<nodes, gps, 2, edim, etype, Physics>: ActionOperator, Physics {

	Integration(size_t interval) { isconst = false; }

	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			double jacobian[4] = { 0, 0, 0, 0 }, inv[4];

			for (size_t n = 0; n < nodes; ++n) {
				jacobian[0] += element.dN[gp][n][0] * element.coords[n][0];
				jacobian[1] += element.dN[gp][n][0] * element.coords[n][1];
				jacobian[2] += element.dN[gp][n][1] * element.coords[n][0];
				jacobian[3] += element.dN[gp][n][1] * element.coords[n][1];
			}

			element.det[gp] = jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2];
			double detJx = 1 / element.det[gp];
			inv[0] =   detJx * jacobian[3];
			inv[1] = - detJx * jacobian[1];
			inv[2] = - detJx * jacobian[2];
			inv[3] =   detJx * jacobian[0];

			for (size_t n = 0; n < nodes; ++n) {
				element.dND[gp][n][0] = inv[0] * element.dN[gp][n][0] + inv[1] * element.dN[gp][n][1];
				element.dND[gp][n][1] = inv[2] * element.dN[gp][n][0] + inv[3] * element.dN[gp][n][1];
			}
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD jacobian0 = zeros(), jacobian1 = zeros(), jacobian2 = zeros(), jacobian3 = zeros();

			for (size_t n = 0; n < nodes; ++n) {
				SIMD coordsX = element.coords[n][0];
				SIMD coordsY = element.coords[n][1];
				SIMD dNX = element.dN[gp][n][0];
				SIMD dNY = element.dN[gp][n][1];

				jacobian0 = jacobian0 + dNX * coordsX;
				jacobian1 = jacobian1 + dNX * coordsY;
				jacobian2 = jacobian2 + dNY * coordsX;
				jacobian3 = jacobian3 + dNY * coordsY;
			}

			SIMD determinant = jacobian0 * jacobian3 - jacobian1 * jacobian2;
			store(element.det[gp], determinant);

			SIMD detJx = ones() / determinant;
			SIMD inv0 =  detJx * jacobian3;
			SIMD inv1 = -detJx * jacobian1;
			SIMD inv2 = -detJx * jacobian2;
			SIMD inv3 =  detJx * jacobian0;

			for (size_t n = 0; n < nodes; ++n) {
				SIMD dNX = load(element.dN[gp][n][0]);
				SIMD dNY = load(element.dN[gp][n][1]);
				element.dND[gp][n][0] = inv0 * dNX + inv1 * dNY;
				element.dND[gp][n][1] = inv2 * dNX + inv3 * dNY;
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct Integration<nodes, gps, 3, edim, etype, Physics>: ActionOperator, Physics {

	Integration(size_t interval) { isconst = false; }

	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			double jacobian[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 }, inv[9];

			for (size_t n = 0; n < nodes; ++n) {
				jacobian[0] += element.dN[gp][n][0] * element.coords[n][0];
				jacobian[1] += element.dN[gp][n][0] * element.coords[n][1];
				jacobian[2] += element.dN[gp][n][0] * element.coords[n][2];
				jacobian[3] += element.dN[gp][n][1] * element.coords[n][0];
				jacobian[4] += element.dN[gp][n][1] * element.coords[n][1];
				jacobian[5] += element.dN[gp][n][1] * element.coords[n][2];
				jacobian[6] += element.dN[gp][n][2] * element.coords[n][0];
				jacobian[7] += element.dN[gp][n][2] * element.coords[n][1];
				jacobian[8] += element.dN[gp][n][2] * element.coords[n][2];
			}
			element.det[gp] =
					+ jacobian[0] * jacobian[4] * jacobian[8]
					+ jacobian[1] * jacobian[5] * jacobian[6]
					+ jacobian[2] * jacobian[3] * jacobian[7]
					- jacobian[2] * jacobian[4] * jacobian[6]
					- jacobian[1] * jacobian[3] * jacobian[8]
					- jacobian[0] * jacobian[5] * jacobian[7];

			double detJx = 1 / element.det[gp];
			inv[0] = detJx * ( jacobian[8] * jacobian[4] - jacobian[7] * jacobian[5]);
			inv[1] = detJx * (-jacobian[8] * jacobian[1] + jacobian[7] * jacobian[2]);
			inv[2] = detJx * ( jacobian[5] * jacobian[1] - jacobian[4] * jacobian[2]);
			inv[3] = detJx * (-jacobian[8] * jacobian[3] + jacobian[6] * jacobian[5]);
			inv[4] = detJx * ( jacobian[8] * jacobian[0] - jacobian[6] * jacobian[2]);
			inv[5] = detJx * (-jacobian[5] * jacobian[0] + jacobian[3] * jacobian[2]);
			inv[6] = detJx * ( jacobian[7] * jacobian[3] - jacobian[6] * jacobian[4]);
			inv[7] = detJx * (-jacobian[7] * jacobian[0] + jacobian[6] * jacobian[1]);
			inv[8] = detJx * ( jacobian[4] * jacobian[0] - jacobian[3] * jacobian[1]);

			for (size_t n = 0; n < nodes; ++n) {
				element.dND[gp][n][0] = inv[0] * element.dN[gp][n][0] + inv[1] * element.dN[gp][n][1] + inv[2] * element.dN[gp][n][2];
				element.dND[gp][n][1] = inv[3] * element.dN[gp][n][0] + inv[4] * element.dN[gp][n][1] + inv[5] * element.dN[gp][n][2];
				element.dND[gp][n][2] = inv[6] * element.dN[gp][n][0] + inv[7] * element.dN[gp][n][1] + inv[8] * element.dN[gp][n][2];
			}
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD jacobian0 = zeros(), jacobian1 = zeros(), jacobian2 = zeros(), jacobian3 = zeros(), jacobian4 = zeros();
			SIMD jacobian5 = zeros(), jacobian6 = zeros(), jacobian7 = zeros(), jacobian8 = zeros();

			for (size_t n = 0; n < nodes; ++n) {
				SIMD coordsX = element.coords[n][0];
				SIMD coordsY = element.coords[n][1];
				SIMD coordsZ = element.coords[n][2];
				SIMD dNX = element.dN[gp][n][0];
				SIMD dNY = element.dN[gp][n][1];
				SIMD dNZ = element.dN[gp][n][2];

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

			SIMD determinant =
					+ jacobian0 * jacobian4 * jacobian8
					+ jacobian1 * jacobian5 * jacobian6
					+ jacobian2 * jacobian3 * jacobian7
					- jacobian2 * jacobian4 * jacobian6
					- jacobian1 * jacobian3 * jacobian8
					- jacobian0 * jacobian5 * jacobian7;
			store(element.det[gp], determinant);

			SIMD detJx = ones() / determinant;
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
				SIMD dNX = element.dN[gp][n][0];
				SIMD dNY = element.dN[gp][n][1];
				SIMD dNZ = element.dN[gp][n][2];
				element.dND[gp][n][0] = inv0 * dNX + inv1 * dNY + inv2 * dNZ;
				element.dND[gp][n][1] = inv3 * dNX + inv4 * dNY + inv5 * dNZ;
				element.dND[gp][n][2] = inv6 * dNX + inv7 * dNY + inv8 * dNZ;
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct Volume: ActionOperator, Physics {
	double &volume;
	Volume(size_t interval, std::vector<double> &volume): volume(volume[interval]) { isconst = false; }

	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			volume += element.det[gp] * element.w[gp];
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD v = element.det[gp] * element.w[gp];
			for (size_t s = 0; s < SIMD::size; ++s) {
				volume += v[s];
			}
		}
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD v = element.det[gp] * element.w[gp];
			for (size_t s = 0; s < size; ++s) {
				volume += v[s];
			}
		}
	}
};

struct BoundaryJacobian: public ActionOperator {
	BoundaryJacobian(
			int interval,
			const ParameterData &coordinates,
			const ParameterData &dN,
			ParameterData &jacobian)
	: coords(coordinates, interval),
	  dN(dN, interval, 0),
	  jacobian(jacobian, interval)
	{ }

	InputParameterIterator coords, dN;
	OutputParameterIterator jacobian;

	void operator++()
	{
		++coords;
		++jacobian;
	}

	void move(int n)
	{
		coords += n;
		jacobian += n;
	}
};

template<size_t nodes, size_t gps>
struct BoundaryFaceJacobian: public BoundaryJacobian {
	using BoundaryJacobian::BoundaryJacobian;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double dND[6] = { 0, 0, 0, 0, 0, 0 };
			M2NMN3<nodes>(1, dN.data + 2 * gpindex * nodes, coords.data, dND);
			double x = dND[1] * dND[5] - dND[2] * dND[4];
			double y = dND[2] * dND[3] - dND[0] * dND[5];
			double z = dND[0] * dND[4] - dND[1] * dND[3];
			jacobian[gpindex] = std::sqrt(x * x + y * y + z * z);
		}
	}
};

template<size_t nodes, size_t gps>
struct BoundaryEdge3DJacobian: public BoundaryJacobian {
	using BoundaryJacobian::BoundaryJacobian;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double dND[3] = { 0, 0, 0 };
			M1NMN3<nodes>(1, dN.data + 1 * gpindex * nodes, coords.data, dND);
			jacobian[gpindex] = std::sqrt(dND[0] * dND[0] + dND[1] * dND[1] + dND[2] * dND[2]);
		}
	}
};

template<size_t nodes, size_t gps>
struct BoundaryEdge2DJacobian: public BoundaryJacobian {
	using BoundaryJacobian::BoundaryJacobian;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double dND[2] = { 0, 0 };
			M1NMN2<nodes>(1, dN.data + 1 * gpindex * nodes, coords.data, dND);
			jacobian[gpindex] = std::sqrt(dND[0] * dND[0] + dND[1] * dND[1]);
		}
	}
};

struct BoundaryNormal: public ActionOperator {
	BoundaryNormal(
			int interval,
			const ParameterData &coordinates,
			const ParameterData &dN,
			ParameterData &jacobian,
			ParameterData &normal)
	: coords(coordinates, interval),
	  dN(dN, interval, 0),
	  jacobian(jacobian, interval),
	  normal(normal, interval)
	{ }

	InputParameterIterator coords, dN;
	OutputParameterIterator jacobian, normal;

	void operator++()
	{
		++coords;
		++jacobian; ++normal;
	}

	void move(int n)
	{
		coords += n;
		jacobian += n;
		normal += n;
	}
};

template<size_t nodes, size_t gps>
struct BoundaryFaceNormal: public BoundaryNormal {
	using BoundaryNormal::BoundaryNormal;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double dND[6] = { 0, 0, 0, 0, 0, 0 };
			M2NMN3<nodes>(1, dN.data + 2 * gpindex * nodes, coords.data, dND);
			double x = dND[1] * dND[5] - dND[2] * dND[4];
			double y = dND[2] * dND[3] - dND[0] * dND[5];
			double z = dND[0] * dND[4] - dND[1] * dND[3];
			jacobian[gpindex] = std::sqrt(x * x + y * y + z * z);
			normal[3 * gpindex + 0] = x / jacobian[gpindex];
			normal[3 * gpindex + 1] = y / jacobian[gpindex];
			normal[3 * gpindex + 2] = z / jacobian[gpindex];
		}
	}
};

template<size_t nodes, size_t gps>
struct BoundaryEdge2DNormal: public BoundaryNormal {
	using BoundaryNormal::BoundaryNormal;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double dND[2] = { 0, 0 };
			M1NMN2<nodes>(1, dN.data + 1 * gpindex * nodes, coords.data, dND);
			jacobian[gpindex] = std::sqrt(dND[0] * dND[0] + dND[1] * dND[1]);
			normal[2 * gpindex + 0] = -dND[1] / jacobian[gpindex];
			normal[2 * gpindex + 1] =  dND[0] / jacobian[gpindex];
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_ */
