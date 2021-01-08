
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_JACOBIAN_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_JACOBIAN_H_

#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"
#include "physics/assembler/math.hpp"

namespace espreso {

struct ElementJacobian: public Operator {
	ElementJacobian(
			const ParameterData &coordinates,
			const ParameterData &dN,
			ParameterData &inversion,
			ParameterData &det,
			ParameterData &dND,
			int interval)
	: Operator(interval, false, Link(interval).inputs(coordinates, dN).outputs(inversion, det, dND)),
	  coords(coordinates, interval, ndim * enodes),
	  dN(dN, interval, 0),
	  inv(inversion, interval, ndim * ndim * egps),
	  det(det, interval, egps),
	  dND(dND, interval, edim * enodes * egps)
	{ }

	InputParameterIterator coords, dN;
	OutputParameterIterator inv, det, dND;

	void operator++()
	{
		++coords;
		++inv; ++det; ++dND;
	}
};

struct ElementJacobian2D: public ElementJacobian {
	GET_NAME(ElementJacobian2D)
	using ElementJacobian::ElementJacobian;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		double jacobian[4] = { 0, 0, 0, 0 };

		for (int n = 0; n < nodes; ++n) {
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
};

struct ElementJacobian3D: public ElementJacobian {
	GET_NAME(ElementJacobian3D)
	using ElementJacobian::ElementJacobian;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		double jacobian[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

		for (int n = 0; n < nodes; ++n) {
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
};

template <class Kernel>
struct ElementIntegration: public ElementOperatorBuilder {
	GET_NAME(ElementIntegration)

	Kernel &kernel;

	ElementIntegration(Kernel &kernel): kernel(kernel)
	{
		kernel.integration.jacobiInversion.addInputs(kernel.coords.node, kernel.integration.dN);
		kernel.integration.jacobiDeterminant.addInputs(kernel.coords.node, kernel.integration.dN);
		kernel.integration.dND.addInputs(kernel.coords.node, kernel.integration.dN);
	}

	void apply(int interval)
	{
		if (info::mesh->dimension == 2) {
			iterate_elements_gps<Kernel>(ElementJacobian2D(kernel.coords.node, kernel.integration.dN, kernel.integration.jacobiInversion, kernel.integration.jacobiDeterminant, kernel.integration.dND, interval));
		}
		if (info::mesh->dimension == 3) {
			iterate_elements_gps<Kernel>(ElementJacobian3D(kernel.coords.node, kernel.integration.dN, kernel.integration.jacobiInversion, kernel.integration.jacobiDeterminant, kernel.integration.dND, interval));
		}
	}
};

struct BoundaryJacobian: public Operator {
	BoundaryJacobian(
			const ParameterData &coordinates,
			const ParameterData &dN,
			ParameterData &jacobian,
			int interval)
	: Operator(interval, jacobian.isconst[interval], Link(interval).inputs(coordinates, dN).outputs(jacobian)),
	  coords(coordinates, interval, ndim * enodes),
	  dN(dN, interval, 0),
	  jacobian(jacobian, interval, egps)
	{ }

	InputParameterIterator coords, dN;
	OutputParameterIterator jacobian;

	void operator++()
	{
		++coords;
		++jacobian;
	}
};

struct BoundaryFaceJacobian: public BoundaryJacobian {
	GET_NAME(BoundaryFaceJacobian)
	using BoundaryJacobian::BoundaryJacobian;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		double dND[6] = { 0, 0, 0, 0, 0, 0 };
		M2NMN3<nodes>(1, dN.data + 2 * gpindex * nodes, coords.data, dND);
		double x = dND[1] * dND[5] - dND[2] * dND[4];
		double y = dND[2] * dND[3] - dND[0] * dND[5];
		double z = dND[0] * dND[4] - dND[1] * dND[3];
		jacobian.data[gpindex] = std::sqrt(x * x + y * y + z * z);
	}
};

struct BoundaryEdge3DJacobian: public BoundaryJacobian {
	GET_NAME(BoundaryEdge3DJacobian)
	using BoundaryJacobian::BoundaryJacobian;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		double dND[3] = { 0, 0, 0 };
		M1NMN3<nodes>(1, dN.data + 1 * gpindex * nodes, coords.data, dND);
		jacobian.data[gpindex] = std::sqrt(dND[0] * dND[0] + dND[1] * dND[1] + dND[2] * dND[2]);
	}
};

struct BoundaryEdge2DJacobian: public BoundaryJacobian {
	GET_NAME(BoundaryEdge2DJacobian)
	using BoundaryJacobian::BoundaryJacobian;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		double dND[2] = { 0, 0 };
		M1NMN2<nodes>(1, dN.data + 1 * gpindex * nodes, coords.data, dND);
		jacobian.data[gpindex] = std::sqrt(dND[0] * dND[0] + dND[1] * dND[1]);
	}
};

template <class Kernel>
struct BoundaryIntegration: public BoundaryOperatorBuilder {
	GET_NAME(BoundaryIntegration)

	Kernel &kernel;

	BoundaryIntegration(Kernel &kernel): kernel(kernel)
	{
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			kernel.integration.boundary.jacobian.regions[r].addInputs(kernel.coords.boundary.node.regions[r], kernel.integration.boundary.dN.regions[r]);
			kernel.integration.boundary.jacobian.regions[r].isset = true;
		}
	}

	void apply(int region, int interval)
	{
		if (info::mesh->boundaryRegions[region]->dimension == 2) {
			iterate_boundary_gps<Kernel>(BoundaryFaceJacobian(kernel.coords.boundary.node.regions[region], kernel.integration.boundary.dN.regions[region], kernel.integration.boundary.jacobian.regions[region], interval), region);
		}
		if (info::mesh->boundaryRegions[region]->dimension == 1) {
			if (info::mesh->dimension == 3) {
				iterate_boundary_gps<Kernel>(BoundaryEdge3DJacobian(kernel.coords.boundary.node.regions[region], kernel.integration.boundary.dN.regions[region], kernel.integration.boundary.jacobian.regions[region], interval), region);
			}
			if (info::mesh->dimension == 2) {
				iterate_boundary_gps<Kernel>(BoundaryEdge2DJacobian(kernel.coords.boundary.node.regions[region], kernel.integration.boundary.dN.regions[region], kernel.integration.boundary.jacobian.regions[region], interval), region);
			}
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_JACOBIAN_H_ */
