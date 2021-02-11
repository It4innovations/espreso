
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_MASS_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_MASS_H_

#include "esinfo/meshinfo.h"
#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

struct MaterialMass: public Operator {
	MaterialMass(ParameterData &density, ParameterData &heatCapacity, ParameterData &thickness, ParameterData &mass, int interval)
	: Operator(interval, false, mass.update[interval]),
	  density(density, interval),
	  heatCapacity(heatCapacity, interval),
	  thickness(thickness, interval),
	  mass(mass, interval)
	{

	}

	InputParameterIterator density, heatCapacity, thickness;
	OutputParameterIterator mass;
};

struct MaterialMass2D: public MaterialMass {
	using MaterialMass::MaterialMass;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		mass[gpindex]  = density[gpindex] * heatCapacity[gpindex] * thickness[gpindex];
	}

	void operator++()
	{
		++density; ++heatCapacity; ++thickness;
		++mass;
	}
};

struct MaterialMass3D: public MaterialMass {
	using MaterialMass::MaterialMass;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		mass[gpindex]  = density[gpindex] * heatCapacity[gpindex];
	}

	void operator++()
	{
		++density; ++heatCapacity;
		++mass;
	}
};

struct MaterialMassBuilder: public ElementOperatorBuilder {
	HeatTransferModuleOpt &kernel;

	MaterialMassBuilder(HeatTransferModuleOpt &kernel): ElementOperatorBuilder("MATERIAL MASS"), kernel(kernel)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		kernel.material.mass.addInput(kernel.material.density);
		kernel.material.mass.addInput(kernel.material.heatCapacity);
		if (info::mesh->dimension == 2) {
			kernel.material.mass.addInput(kernel.thickness.gp);
		}
		kernel.material.mass.resize();
		kernel.addParameter(kernel.material.mass);
		return true;
	}

	void apply(int interval)
	{
		if (info::mesh->dimension == 2) {
			iterate_elements_gps<HeatTransferModuleOpt::NGP>(MaterialMass2D(kernel.material.density, kernel.material.heatCapacity, kernel.thickness.gp, kernel.material.mass, interval));
		}
		if (info::mesh->dimension == 3) {
			iterate_elements_gps<HeatTransferModuleOpt::NGP>(MaterialMass3D(kernel.material.density, kernel.material.heatCapacity, kernel.thickness.gp, kernel.material.mass, interval));
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_MASS_H_ */
