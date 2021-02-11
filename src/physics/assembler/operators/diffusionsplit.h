
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_DIFFUSIONSPLIT_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_DIFFUSIONSPLIT_H_

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

struct DiffusionSplitBase: public Operator {
	DiffusionSplitBase(ElementData *gradient, ParameterData &dND, ParameterData &conductivity, ParameterData &mass, ParameterData &xi, int interval)
	: Operator(interval, false, xi.update[interval]),
	  gradient(gradient->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin, info::mesh->dimension),
	  dND(dND, interval),
	  conductivity(conductivity, interval), // we touch only the first element in both isotropic and non-isotropic cases
	  mass(mass, interval),
	  xi(xi, interval)
	{

	}

	const double C1 = 1, C2 = 6;
	InputParameterIterator gradient, dND, conductivity, mass;
	OutputParameterIterator xi;

	void operator++()
	{
		++gradient; ++dND; ++conductivity; ++mass;
		++xi;
	}
};

struct DiffusionSplit2D: public DiffusionSplitBase {
	using DiffusionSplitBase::DiffusionSplitBase;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		double g[nodes];
		M12M2N<nodes>(gradient.data, dND.data + 2 * gpindex * nodes, g);

		double grad = std::sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]);
		double norm = 0;
		for (int n = 0; n < nodes; ++n) {
			norm += g[n] * g[n];
		}
		double H = 2 * std::sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]) / norm;
		double T = (C1 * H * H) / (conductivity[0] * C2 + H * H * (mass[0] / step::time::shift));
		xi[0] = std::max(1., 1 / (1 - T * mass[0] / step::time::shift));
	}
};

struct DiffusionSplit3D: public DiffusionSplitBase {
	using DiffusionSplitBase::DiffusionSplitBase;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		double g[nodes];
		M12M2N<nodes>(gradient.data, dND.data + 2 * gpindex * nodes, g);

		double grad = std::sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]);
		double norm = 0;
		for (int n = 0; n < nodes; ++n) {
			norm += g[n] * g[n];
		}
		double H = 2 * std::sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]) / norm;
		double T = (C1 * H * H) / (conductivity[0] * C2 + H * H * (mass[0] / step::time::shift));
		xi[0] = std::max(1., 1 / (1 - T * mass[0] / step::time::shift));
	}
};

struct DiffusionSplit: public ElementOperatorBuilder {
	HeatTransferModuleOpt &kernel;

	DiffusionSplit(HeatTransferModuleOpt &kernel): ElementOperatorBuilder("DIFFUSION SPLIT"), kernel(kernel)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		kernel.gradient.xi.addInput(kernel.gradient.output);
		kernel.gradient.xi.addInput(kernel.integration.dND);
		kernel.gradient.xi.addInput(kernel.material.conductivity);
		kernel.gradient.xi.addInput(kernel.material.mass);
		kernel.gradient.xi.resize();
		kernel.addParameter(kernel.gradient.xi);
		return true;
	}

	void apply(int interval)
	{
//		if (info::mesh->dimension == 2) {
//			iterate_elements_gps<HeatTransferModuleOpt>(HeatGradient2D(kernel.integration.dND, kernel.temp.node, kernel.gradient.output, interval));
//		}
//		if (info::mesh->dimension == 3) {
//			iterate_elements_gps<HeatTransferModuleOpt>(HeatGradient3D(kernel.integration.dND, kernel.temp.node, kernel.gradient.output, interval));
//		}
	}
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_DIFFUSIONSPLIT_H_ */
