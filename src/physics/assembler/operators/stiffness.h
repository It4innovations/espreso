
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
			const ParameterData &thickness,
			ParameterData &stiffness,
			int interval)
	: Operator(interval, stiffness.isconst[interval], Link(interval).inputs(dND, weight, determinant, conductivity, thickness).outputs(stiffness)),
	  dND(dND, interval, dND.size),
	  weight(weight, interval, 0),
	  determinant(determinant, interval, egps),
	  conductivity(conductivity, interval, conductivity.size),
	  thickness(thickness, interval, 1),
	  stiffness(stiffness, interval, stiffness.size)
	{
		if (update) {
			std::fill((stiffness.data->begin() + interval)->data(), (stiffness.data->begin() + interval + 1)->data(), 0);
		}
	}

	InputParameterIterator dND, weight, determinant, conductivity, thickness;
	OutputParameterIterator stiffness;

	void operator++()
	{
		++dND; ++determinant; ++conductivity; ++thickness;
		++stiffness;
	}
};

struct Stiffness2DHeatIsotropic: public Stiffness {
	GET_NAME(Stiffness2DHeatIsotropic)
	using Stiffness::Stiffness;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		ADDMN2M2N<nodes>(thickness[gpindex] * determinant[gpindex] * weight[gpindex] * conductivity[gpindex], dND.data + 2 * nodes * gpindex, stiffness.data);
	}
};

struct Stiffness2DHeat: public Stiffness {
	GET_NAME(Stiffness2DHeat)
	using Stiffness::Stiffness;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		ADDMN2M22M2N<nodes>(thickness[gpindex] * determinant[gpindex] * weight[gpindex], conductivity.data + 4 * gpindex, dND.data + 2 * nodes * gpindex, stiffness.data);
	}
};

struct Stiffness3DHeatIsotropic: public Stiffness {
	GET_NAME(Stiffness3DHeatIsotropic)
	using Stiffness::Stiffness;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		ADDMN3M3N<nodes>(determinant[gpindex] * weight[gpindex] * conductivity[gpindex], dND.data + 3 * nodes * gpindex, stiffness.data);
	}
};

struct Stiffness3DHeat: public Stiffness {
	GET_NAME(Stiffness3DHeat)
	using Stiffness::Stiffness;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		KMN3M33M3N<nodes>(determinant[gpindex] * weight[gpindex], conductivity.data + 9 * gpindex, dND.data + 3 * nodes * gpindex, stiffness.data);
	}
};

struct HeatStiffness: public ElementOperatorBuilder {
	GET_NAME(HeatStiffness)
	HeatTransferModuleOpt &kernel;

	HeatStiffness(HeatTransferModuleOpt &kernel): kernel(kernel)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		if (info::mesh->dimension == 2) {
			kernel.linearSystem.stiffness.addInputs(kernel.integration.dND, kernel.integration.weight, kernel.integration.jacobiDeterminant, kernel.material.conductivity, kernel.thickness.gp);
		}
		if (info::mesh->dimension == 3) {
			kernel.linearSystem.stiffness.addInputs(kernel.integration.dND, kernel.integration.weight, kernel.integration.jacobiDeterminant, kernel.material.conductivity);
		}
		return true;
	}

	void apply(int interval)
	{
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
		if (mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
			if (info::mesh->dimension == 2) {
				iterate_elements_gps<HeatTransferModuleOpt>(Stiffness2DHeatIsotropic(kernel.integration.dND, kernel.integration.weight, kernel.integration.jacobiDeterminant, kernel.material.conductivityIsotropic, kernel.thickness.gp, kernel.linearSystem.stiffness, interval));
			}
			if (info::mesh->dimension == 3) {
				iterate_elements_gps<HeatTransferModuleOpt>(Stiffness3DHeatIsotropic(kernel.integration.dND, kernel.integration.weight, kernel.integration.jacobiDeterminant, kernel.material.conductivityIsotropic, kernel.thickness.gp, kernel.linearSystem.stiffness, interval));
			}
		} else {
			if (info::mesh->dimension == 2) {
				iterate_elements_gps<HeatTransferModuleOpt>(Stiffness2DHeat(kernel.integration.dND, kernel.integration.weight, kernel.integration.jacobiDeterminant, kernel.material.conductivity, kernel.linearSystem.stiffness, kernel.thickness.gp, interval));
			}
			if (info::mesh->dimension == 3) {
				iterate_elements_gps<HeatTransferModuleOpt>(Stiffness3DHeat(kernel.integration.dND, kernel.integration.weight, kernel.integration.jacobiDeterminant, kernel.material.conductivity, kernel.linearSystem.stiffness, kernel.thickness.gp, interval));
			}
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_STIFFNESS_H_ */
