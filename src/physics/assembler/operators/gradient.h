
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_GRADIENT_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_GRADIENT_H_

#include "mesh/store/elementstore.h"
#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

struct OutputGradient: public Operator {
	OutputGradient(ParameterData &dND, ParameterData &temperature, ElementData *gradient, int interval)
	: Operator(interval, false, true),
	  dND(dND, interval),
	  temp(temperature, interval),
	  gradient(gradient->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin, info::mesh->dimension)
	{
		if (update) {
			std::fill(gradient->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin, gradient->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].end, 0);
		}
	}

	InputParameterIterator dND, temp;
	OutputParameterIterator gradient;

	void operator++()
	{
		++dND; ++temp;
		++gradient;
	}
};

struct HeatGradient2D: public OutputGradient {
	using OutputGradient::OutputGradient;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		ADDM2NMN1<nodes>(1. / gps, dND.data + 2 * nodes * gpindex, temp.data, gradient.data);
	}
};

struct HeatGradient3D: public OutputGradient {
	using OutputGradient::OutputGradient;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		ADDM3NMN1<nodes>(1. / gps, dND.data + 3 * nodes * gpindex, temp.data, gradient.data);
	}
};

struct Gradient: public ElementOperatorBuilder {
	HeatTransferModuleOpt &kernel;

	Gradient(HeatTransferModuleOpt &kernel): ElementOperatorBuilder("ELEMENTS GRADIENT"), kernel(kernel)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		return true;
	}

	void apply(int interval)
	{
		if (info::mesh->dimension == 2) {
			iterate_elements_gps<HeatTransferModuleOpt::NGP>(HeatGradient2D(kernel.integration.dND, kernel.temp.node, kernel.gradient.output, interval));
		}
		if (info::mesh->dimension == 3) {
			iterate_elements_gps<HeatTransferModuleOpt::NGP>(HeatGradient3D(kernel.integration.dND, kernel.temp.node, kernel.gradient.output, interval));
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_GRADIENT_H_ */
