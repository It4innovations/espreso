
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_RHS_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_RHS_H_

#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"
#include "physics/assembler/math.hpp"

namespace espreso {

//struct HeatFlow: public Operator {
//	GET_NAME(HeatFlow)
//
//	HeatFlow(double area, const ParameterData &heatFlow, ParameterData &q, int interval)
//	: Operator(interval, q.isconst[interval], Link(interval).inputs(heatFlow).outputs(q)),
//	  area(area), heatFlow(heatFlow, interval, egps), q(q, interval, egps) { }
//
//	double area;
//	InputParameterIterator heatFlow;
//	OutputParameterIterator q;
//
//	template<int nodes, int gps>
//	void operator()(int gpindex)
//	{
//		q.data[gpindex] += heatFlow.data[gpindex] / area;
//	}
//
//	void operator++()
//	{
//		++heatFlow;
//	}
//};
//
//struct HeatFlux: public Operator {
//	GET_NAME(HeatFlux)
//
//	HeatFlux(const ParameterData &heatFlux, ParameterData &q, int interval)
//	: Operator(interval, q.isconst[interval], Link(interval).inputs(heatFlux).outputs(q)),
//	  heatFlux(heatFlux, interval, egps), q(q, interval, egps) { }
//
//	InputParameterIterator heatFlux;
//	OutputParameterIterator q;
//
//	template<int nodes, int gps>
//	void operator()(int gpindex)
//	{
//		q.data[gpindex] += heatFlux.data[gpindex];
//	}
//
//	void operator++()
//	{
//		++heatFlux;
//	}
//};
//
//struct ConvectionInitial: public Operator {
//	GET_NAME(ConvectionInitial)
//
//	ConvectionInitial(const ParameterData &htc, const ParameterData &extTemp, ParameterData &q, int interval)
//	: Operator(interval, q.isconst[interval], Link(interval).inputs(htc, extTemp).outputs(q)),
//	  htc(htc, interval, egps), extTemp(extTemp, interval, egps), q(q, interval, egps) { }
//
//	InputParameterIterator htc, extTemp;
//	OutputParameterIterator q;
//
//	template<int nodes, int gps>
//	void operator()(int gpindex)
//	{
//		q.data[gpindex] += htc.data[gpindex] * extTemp.data[gpindex];
//	}
//
//	void operator++()
//	{
//		++htc; ++extTemp;
//		++q;
//	}
//};

struct HeatQ: public Operator {
	HeatQ(double area, const ParameterData &heatFlow, const ParameterData &heatFlux, const ParameterData &htc, const ParameterData &extTemp, ParameterData &q, int interval)
	: Operator(interval, q.isconst[interval], q.update[interval]),
	  area(area), heatFlow(heatFlow, interval),
	  heatFlux(heatFlux, interval),
	  htc(htc, interval), extTemp(extTemp, interval),
	  q(q, interval)
	{
		if (update) {
			std::fill((q.data->begin() + interval)->data(), (q.data->begin() + interval + 1)->data(), 0);
		}
	}

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		q.data[gpindex] += heatFlow.data[gpindex] / area;
		q.data[gpindex] += heatFlux.data[gpindex];
		q.data[gpindex] += htc.data[gpindex] * extTemp.data[gpindex];
	}

	void operator++()
	{
		++heatFlow; ++heatFlux;
		++htc; ++extTemp;
		++q;
	}

	double area;
	InputParameterIterator heatFlow, heatFlux, htc, extTemp;
	OutputParameterIterator q;
};

struct HeatRHS2D: public Operator {
	HeatRHS2D(const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &thickness, const ParameterData &q, ParameterData &rhs, int interval)
	: Operator(interval, rhs.isconst[interval], rhs.update[interval]),
	  N(N, interval),
	  weight(weight, interval),
	  J(J, interval),
	  thickness(thickness, interval),
	  q(q, interval),
	  rhs(rhs, interval)
	{ }

	InputParameterIterator N, weight, J;
	InputParameterIterator thickness, q;
	OutputParameterIterator rhs;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		for (int n = 0; n < nodes; ++n) {
			rhs.data[n] += thickness.data[gpindex] * J.data[gpindex] * weight.data[gpindex] * q.data[gpindex] * N.data[gpindex * nodes + n];
		}
	}

	void operator++()
	{
		++weight; ++J;
		++thickness; ++q;
		++rhs;
	}
};

struct HeatRHS3D: public Operator {
	HeatRHS3D(const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &q, ParameterData &rhs, int interval)
	: Operator(interval, rhs.isconst[interval], rhs.update[interval]),
	  N(N, interval),
	  weight(weight, interval),
	  J(J, interval),
	  q(q, interval),
	  rhs(rhs, interval)
	{
		if (update) {
			std::fill((rhs.data->begin() + interval)->data(), (rhs.data->begin() + interval + 1)->data(), 0);
		}
	}

	InputParameterIterator N, weight, J;
	InputParameterIterator q;
	OutputParameterIterator rhs;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		for (int n = 0; n < nodes; ++n) {
			rhs.data[n] += J.data[gpindex] * weight.data[gpindex] * q.data[gpindex] * N.data[gpindex * nodes + n];
		}
	}

	void operator++()
	{
		++weight; ++J;
		++q;
		++rhs;
	}
};

struct HeatRHS: public BoundaryOperatorBuilder {
	HeatTransferModuleOpt &kernel;

	HeatRHS(HeatTransferModuleOpt &kernel): BoundaryOperatorBuilder("HEAT TRANSFER RHS"), kernel(kernel)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				if (kernel.heatFlow.gp.regions[r].data == NULL) {
					kernel.heatFlow.gp.regions[r].resize();
					kernel.addParameter(kernel.heatFlow.gp.regions[r]);
				}
				if (kernel.heatFlux.gp.regions[r].data == NULL) {
					kernel.heatFlux.gp.regions[r].resize();
					kernel.addParameter(kernel.heatFlux.gp.regions[r]);
				}
				if (kernel.convection.heatTransferCoeficient.gp.regions[r].data == NULL) {
					kernel.convection.heatTransferCoeficient.gp.regions[r].resize();
					kernel.convection.externalTemperature.gp.regions[r].resize();
					kernel.addParameter(kernel.convection.externalTemperature.gp.regions[r]);
				}

				if (info::mesh->dimension == 2) {
					kernel.elements.boundary.rhs.regions[r].addInput(kernel.thickness.boundary.gp.regions[r]);
				}
				kernel.q.gp.regions[r].addInput(kernel.heatFlow.gp.regions[r]);
				kernel.q.gp.regions[r].addInput(kernel.heatFlux.gp.regions[r]);
				kernel.q.gp.regions[r].addInput(kernel.convection.heatTransferCoeficient.gp.regions[r]);
				kernel.q.gp.regions[r].addInput(kernel.convection.externalTemperature.gp.regions[r]);
				kernel.q.gp.regions[r].resize();
				kernel.addParameter(kernel.q.gp.regions[r]);

				kernel.elements.boundary.rhs.regions[r].addInput(kernel.q.gp.regions[r]);
				kernel.elements.boundary.rhs.regions[r].addInput(kernel.integration.boundary.jacobian.regions[r]);
				kernel.elements.boundary.rhs.regions[r].addInput(kernel.integration.boundary.weight.regions[r]);
				kernel.elements.boundary.rhs.regions[r].resize();
				kernel.addParameter(kernel.elements.boundary.rhs.regions[r]);
			}
		}
		return true;
	}

	void apply(int region, int interval)
	{
		iterate_boundary_gps<HeatTransferModuleOpt::NGP>(HeatQ(
				info::mesh->boundaryRegions[region]->area, kernel.heatFlow.gp.regions[region],
				kernel.heatFlux.gp.regions[region],
				kernel.convection.heatTransferCoeficient.gp.regions[region], kernel.convection.externalTemperature.gp.regions[region],
				kernel.q.gp.regions[region], interval), region);

		if (info::mesh->dimension == 2) {
			iterate_boundary_gps<HeatTransferModuleOpt::NGP>(HeatRHS2D(
					kernel.integration.boundary.N.regions[region], kernel.integration.boundary.weight.regions[region], kernel.integration.boundary.jacobian.regions[region],
					kernel.thickness.boundary.gp.regions[region], kernel.q.gp.regions[region],
					kernel.elements.boundary.rhs.regions[region], interval), region);
		}
		if (info::mesh->dimension == 3) {
			iterate_boundary_gps<HeatTransferModuleOpt::NGP>(HeatRHS3D(
					kernel.integration.boundary.N.regions[region], kernel.integration.boundary.weight.regions[region], kernel.integration.boundary.jacobian.regions[region],
					kernel.q.gp.regions[region], kernel.elements.boundary.rhs.regions[region], interval), region);
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_RHS_H_ */
