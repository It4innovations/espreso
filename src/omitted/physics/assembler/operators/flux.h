
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_FLUX_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_FLUX_H_

#include "mesh/store/elementstore.h"
#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

struct OutputFlux: public Operator {
	OutputFlux(ParameterData &dND, ParameterData &temperature, ParameterData &conductivity, ElementData *flux, int interval)
	: Operator(interval, false, true),
	  dND(dND, interval),
	  temp(temperature, interval),
	  conductivity(conductivity, interval),
	  flux(flux->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin, info::mesh->dimension)
	{
		if (update) {
			std::fill(flux->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin, flux->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].end, 0);
		}
	}

	InputParameterIterator dND, temp, conductivity;
	OutputParameterIterator flux;

	void operator++()
	{
		++dND; ++temp; ++conductivity;
		++flux;
	}
};

struct OutputFluxIsotropic2D: public OutputFlux {
	using OutputFlux::OutputFlux;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		ADDM2NMN1<nodes>(conductivity.data[gpindex] / gps, dND.data + 2 * nodes * gpindex, temp.data, flux.data);
	}
};

struct OutputFluxIsotropic3D: public OutputFlux {
	using OutputFlux::OutputFlux;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		ADDM2NMN1<nodes>(conductivity.data[gpindex] / gps, dND.data + 3 * nodes * gpindex, temp.data, flux.data);
	}
};

struct OutputFlux2D: public OutputFlux {
	using OutputFlux::OutputFlux;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		ADDM22M2NMN1<nodes>(1. / gps, conductivity.data + 4 * gpindex, dND.data + 2 * nodes * gpindex, temp.data, flux.data);
	}
};

struct OutputFlux3D: public OutputFlux {
	using OutputFlux::OutputFlux;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		ADDM33M3NMN1<nodes>(1. / gps, conductivity.data + 9 * gpindex, dND.data + 3 * nodes * gpindex, temp.data, flux.data);
	}
};

struct Flux: public ElementOperatorBuilder {
	HeatTransferModuleOpt &kernel;

	Flux(HeatTransferModuleOpt &kernel): ElementOperatorBuilder("ELEMENTS FLUX"), kernel(kernel)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		return true;
	}

	void apply(int interval)
	{
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
		if (info::mesh->dimension == 2) {
			if (mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
				iterate_elements_gps<HeatTransferModuleOpt::NGP>(OutputFluxIsotropic2D(kernel.integration.dND, kernel.temp.node, kernel.material.conductivityIsotropic, kernel.flux.output, interval));
			} else {
				iterate_elements_gps<HeatTransferModuleOpt::NGP>(OutputFlux2D(kernel.integration.dND, kernel.temp.node, kernel.material.conductivity, kernel.flux.output, interval));
			}
		}
		if (info::mesh->dimension == 3) {
			if (mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
				iterate_elements_gps<HeatTransferModuleOpt::NGP>(OutputFluxIsotropic3D(kernel.integration.dND, kernel.temp.node, kernel.material.conductivityIsotropic, kernel.flux.output, interval));
			} else {
				iterate_elements_gps<HeatTransferModuleOpt::NGP>(OutputFlux3D(kernel.integration.dND, kernel.temp.node, kernel.material.conductivity, kernel.flux.output, interval));
			}
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_FLUX_H_ */
