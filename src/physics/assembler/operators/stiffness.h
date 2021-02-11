
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_STIFFNESS_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_STIFFNESS_H_

#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"
#include "physics/assembler/math.hpp"

namespace espreso {

struct Stiffness: public Operator {
	Stiffness(
			const ParameterData &dND,
			const ParameterData &weight,
			const ParameterData &determinant,
			const ParameterData &conductivity,
			const ParameterData &xi,
			const ParameterData &thickness,
			ParameterData &stiffness,
			int interval)
	: Operator(interval, stiffness.isconst[interval], stiffness.update[interval]),
	  dND(dND, interval),
	  weight(weight, interval, 0),
	  determinant(determinant, interval),
	  conductivity(conductivity, interval),
	  xi(xi, interval),
	  thickness(thickness, interval),
	  stiffness(stiffness, interval)
	{
		if (update) {
			std::fill((stiffness.data->begin() + interval)->data(), (stiffness.data->begin() + interval + 1)->data(), 0);
		}
	}

	InputParameterIterator dND, weight, determinant, conductivity, xi, thickness;
	OutputParameterIterator stiffness;

	void operator++()
	{
		++dND; ++determinant; ++conductivity; ++xi; ++thickness;
		++stiffness;
	}
};

struct Stiffness2DHeatIsotropic: public Stiffness {
	using Stiffness::Stiffness;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		ADDMN2M2N<nodes>(thickness[gpindex] * xi[gpindex] * determinant[gpindex] * weight[gpindex] * conductivity[gpindex], dND.data + 2 * nodes * gpindex, stiffness.data);
	}
};

struct Stiffness2DHeat: public Stiffness {
	using Stiffness::Stiffness;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		ADDMN2M22M2N<nodes>(thickness[gpindex] * xi[gpindex] * determinant[gpindex] * weight[gpindex], conductivity.data + 4 * gpindex, dND.data + 2 * nodes * gpindex, stiffness.data);
	}
};

struct Stiffness3DHeatIsotropic: public Stiffness {
	using Stiffness::Stiffness;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		ADDMN3M3N<nodes>(xi[gpindex] * determinant[gpindex] * weight[gpindex] * conductivity[gpindex], dND.data + 3 * nodes * gpindex, stiffness.data);
	}
};

struct Stiffness3DHeat: public Stiffness {
	using Stiffness::Stiffness;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		KMN3M33M3N<nodes>(xi[gpindex] * determinant[gpindex] * weight[gpindex], conductivity.data + 9 * gpindex, dND.data + 3 * nodes * gpindex, stiffness.data);
	}
};

struct HeatStiffness: public ElementOperatorBuilder {
	HeatTransferModuleOpt &kernel;

	HeatStiffness(HeatTransferModuleOpt &kernel): ElementOperatorBuilder("HEAT TRANSFER STIFFNESS"), kernel(kernel)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		if (info::mesh->dimension == 2) {
			kernel.elements.stiffness.addInput(kernel.thickness.gp);
		}
		kernel.elements.stiffness.addInput(kernel.integration.dND);
		kernel.elements.stiffness.addInput(kernel.integration.weight);
		kernel.elements.stiffness.addInput(kernel.integration.jacobiDeterminant);
		kernel.elements.stiffness.addInput(kernel.material.conductivity);
		kernel.elements.stiffness.addInput(kernel.gradient.xi);
		kernel.elements.stiffness.resize();
		kernel.addParameter(kernel.elements.stiffness);
		return true;
	}

	void apply(int interval)
	{
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
		if (mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
			if (info::mesh->dimension == 2) {
				iterate_elements_gps<HeatTransferModuleOpt::NGP>(Stiffness2DHeatIsotropic(kernel.integration.dND, kernel.integration.weight, kernel.integration.jacobiDeterminant, kernel.material.conductivityIsotropic, kernel.gradient.xi, kernel.thickness.gp, kernel.elements.stiffness, interval));
			}
			if (info::mesh->dimension == 3) {
				iterate_elements_gps<HeatTransferModuleOpt::NGP>(Stiffness3DHeatIsotropic(kernel.integration.dND, kernel.integration.weight, kernel.integration.jacobiDeterminant, kernel.material.conductivityIsotropic, kernel.gradient.xi, kernel.thickness.gp, kernel.elements.stiffness, interval));
			}
		} else {
			if (info::mesh->dimension == 2) {
				iterate_elements_gps<HeatTransferModuleOpt::NGP>(Stiffness2DHeat(kernel.integration.dND, kernel.integration.weight, kernel.integration.jacobiDeterminant, kernel.material.conductivity, kernel.gradient.xi, kernel.thickness.gp, kernel.elements.stiffness, interval));
			}
			if (info::mesh->dimension == 3) {
				iterate_elements_gps<HeatTransferModuleOpt::NGP>(Stiffness3DHeat(kernel.integration.dND, kernel.integration.weight, kernel.integration.jacobiDeterminant, kernel.material.conductivity, kernel.gradient.xi, kernel.thickness.gp, kernel.elements.stiffness, interval));
			}
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_STIFFNESS_H_ */
