
#include "integration.h"

#include "analysis/assembler/math.hpp"
#include "esinfo/meshinfo.h"

using namespace espreso;

struct BoundaryJacobian: public Operator {
	BoundaryJacobian(
			const ParameterData &coordinates,
			const ParameterData &dN,
			ParameterData &jacobian,
			int interval)
	: Operator(interval, jacobian.isconst[interval], jacobian.update[interval]),
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
};

struct BoundaryFaceJacobian: public BoundaryJacobian {
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
	using BoundaryJacobian::BoundaryJacobian;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		double dND[2] = { 0, 0 };
		M1NMN2<nodes>(1, dN.data + 1 * gpindex * nodes, coords.data, dND);
		jacobian.data[gpindex] = std::sqrt(dND[0] * dND[0] + dND[1] * dND[1]);
	}
};


//BoundaryIntegration::BoundaryIntegration(AX_HeatTransfer &module): BoundaryOperatorBuilder("BOUNDARY INTERGRATION")
//{
//	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
//		module.integration.boundary.jacobian.regions[r].addInput(module.coords.boundary.node.regions[r]);
//		module.integration.boundary.jacobian.regions[r].addInput(module.integration.boundary.dN.regions[r]);
//		module.integration.boundary.jacobian.regions[r].isset = true;
//		module.integration.boundary.jacobian.regions[r].resize();
//		module.addParameter(module.integration.boundary.jacobian.regions[r]);
//	}
//}
//
//void BoundaryIntegration::apply(int region, int interval)
//{
////	if (info::mesh->boundaryRegions[region]->dimension == 2) {
////		iterate_boundary_gps<typename Assembler::NGP>(BoundaryFaceJacobian(module.coords.boundary.node.regions[region], module.integration.boundary.dN.regions[region], module.integration.boundary.jacobian.regions[region], interval), region);
////	}
////	if (info::mesh->boundaryRegions[region]->dimension == 1) {
////		if (info::mesh->dimension == 3) {
////			iterate_boundary_gps<typename Assembler::NGP>(BoundaryEdge3DJacobian(module.coords.boundary.node.regions[region], module.integration.boundary.dN.regions[region], module.integration.boundary.jacobian.regions[region], interval), region);
////		}
////		if (info::mesh->dimension == 2) {
////			iterate_boundary_gps<typename Assembler::NGP>(BoundaryEdge2DJacobian(module.coords.boundary.node.regions[region], module.integration.boundary.dN.regions[region], module.integration.boundary.jacobian.regions[region], interval), region);
////		}
////	}
//}
