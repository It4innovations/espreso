
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

namespace espreso {

struct ElementJacobian: public ActionOperator {
	ElementJacobian(
			int interval,
			const ParameterData &coordinates,
			const ParameterData &dN,
			ParameterData &inversion,
			ParameterData &det,
			ParameterData &dND)
	: ActionOperator(interval, false, inversion.update[interval] || det.update[interval] || dND.update[interval]),
	  coords(coordinates, interval),
	  dN(dN, interval, 0),
	  inv(inversion, interval),
	  det(det, interval),
	  dND(dND, interval)
	{

	}

	InputParameterIterator coords, dN;
	OutputParameterIterator inv, det, dND;

	void operator++()
	{
		++coords;
		++inv; ++det; ++dND;
	}

	ElementJacobian& operator+=(const size_t rhs)
	{
		coords +=rhs;
		inv += rhs; det += rhs; dND += rhs;
		return *this;
	}

	void reset()
	{

	}
};

template<size_t nodes, size_t gps>
struct ElementJacobian2D: public ElementJacobian {
	using ElementJacobian::ElementJacobian;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double jacobian[4] = { 0, 0, 0, 0 };

			for (size_t n = 0; n < nodes; ++n) {
				jacobian[0] += dN[2 * gpindex * nodes + n + 0 * nodes] * coords[2 * n + 0];
				jacobian[1] += dN[2 * gpindex * nodes + n + 0 * nodes] * coords[2 * n + 1];
				jacobian[2] += dN[2 * gpindex * nodes + n + 1 * nodes] * coords[2 * n + 0];
				jacobian[3] += dN[2 * gpindex * nodes + n + 1 * nodes] * coords[2 * n + 1];
			}

			det[gpindex] = jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2];
			double detJx = 1 / det[gpindex];
			inv[4 * gpindex + 0] =   detJx * jacobian[3];
			inv[4 * gpindex + 1] = - detJx * jacobian[1];
			inv[4 * gpindex + 2] = - detJx * jacobian[2];
			inv[4 * gpindex + 3] =   detJx * jacobian[0];

			M22M2N<nodes>(inv.data + 4 * gpindex, dN.data + 2 * gpindex * nodes, dND.data + 2 * gpindex * nodes);
		}
	}
};

template<size_t nodes, size_t gps>
struct ElementJacobian3D: public ElementJacobian {
	using ElementJacobian::ElementJacobian;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double jacobian[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

			for (size_t n = 0; n < nodes; ++n) {
				jacobian[0] += dN[3 * gpindex * nodes + n + 0 * nodes] * coords[3 * n + 0];
				jacobian[1] += dN[3 * gpindex * nodes + n + 0 * nodes] * coords[3 * n + 1];
				jacobian[2] += dN[3 * gpindex * nodes + n + 0 * nodes] * coords[3 * n + 2];
				jacobian[3] += dN[3 * gpindex * nodes + n + 1 * nodes] * coords[3 * n + 0];
				jacobian[4] += dN[3 * gpindex * nodes + n + 1 * nodes] * coords[3 * n + 1];
				jacobian[5] += dN[3 * gpindex * nodes + n + 1 * nodes] * coords[3 * n + 2];
				jacobian[6] += dN[3 * gpindex * nodes + n + 2 * nodes] * coords[3 * n + 0];
				jacobian[7] += dN[3 * gpindex * nodes + n + 2 * nodes] * coords[3 * n + 1];
				jacobian[8] += dN[3 * gpindex * nodes + n + 2 * nodes] * coords[3 * n + 2];
			}
			det[gpindex] =
					+ jacobian[0] * jacobian[4] * jacobian[8]
					+ jacobian[1] * jacobian[5] * jacobian[6]
					+ jacobian[2] * jacobian[3] * jacobian[7]
					- jacobian[2] * jacobian[4] * jacobian[6]
					- jacobian[1] * jacobian[3] * jacobian[8]
					- jacobian[0] * jacobian[5] * jacobian[7];

			double detJx = 1 / det[gpindex];
			inv[9 * gpindex + 0] = detJx * ( jacobian[8] * jacobian[4] - jacobian[7] * jacobian[5]);
			inv[9 * gpindex + 1] = detJx * (-jacobian[8] * jacobian[1] + jacobian[7] * jacobian[2]);
			inv[9 * gpindex + 2] = detJx * ( jacobian[5] * jacobian[1] - jacobian[4] * jacobian[2]);
			inv[9 * gpindex + 3] = detJx * (-jacobian[8] * jacobian[3] + jacobian[6] * jacobian[5]);
			inv[9 * gpindex + 4] = detJx * ( jacobian[8] * jacobian[0] - jacobian[6] * jacobian[2]);
			inv[9 * gpindex + 5] = detJx * (-jacobian[5] * jacobian[0] + jacobian[3] * jacobian[2]);
			inv[9 * gpindex + 6] = detJx * ( jacobian[7] * jacobian[3] - jacobian[6] * jacobian[4]);
			inv[9 * gpindex + 7] = detJx * (-jacobian[7] * jacobian[0] + jacobian[6] * jacobian[1]);
			inv[9 * gpindex + 8] = detJx * ( jacobian[4] * jacobian[0] - jacobian[3] * jacobian[1]);

			M33M3N<nodes>(inv.data + 9 * gpindex, dN.data + 3 * gpindex * nodes, dND.data + 3 * gpindex * nodes);
		}
	}
};

struct BoundaryJacobian: public ActionOperator {
	BoundaryJacobian(
			int interval,
			const ParameterData &coordinates,
			const ParameterData &dN,
			ParameterData &jacobian)
	: ActionOperator(interval, jacobian.isconst[interval], jacobian.update[interval]),
	  coords(coordinates, interval),
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

	void reset()
	{

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
			jacobian.data[gpindex] = std::sqrt(x * x + y * y + z * z);
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
			jacobian.data[gpindex] = std::sqrt(dND[0] * dND[0] + dND[1] * dND[1] + dND[2] * dND[2]);
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
			jacobian.data[gpindex] = std::sqrt(dND[0] * dND[0] + dND[1] * dND[1]);
		}
	}
};


}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_ */