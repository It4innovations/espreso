
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_DIFFUSIONSPLIT_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_DIFFUSIONSPLIT_H_

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

struct DiffusionSplitBase: public Operator {
	DiffusionSplitBase(ElementData *gradient, ParameterData &dND, ParameterData &conductivity, ParameterData &mass, ParameterData &xi, int interval)
	: Operator(interval, false, Link(interval).resultIn(gradient).inputs(dND, conductivity, mass).outputs(xi)),
	  gradient(gradient->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin, info::mesh->dimension),
	  dND(dND, interval, edim * enodes * egps),
	  conductivity(conductivity, interval, conductivity.size), // we touch only the first element in both isotropic and non-isotropic cases
	  mass(mass, interval, egps),
	  xi(xi, interval, egps)
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
	GET_NAME(DiffusionSplit2D)
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
	GET_NAME(DiffusionSplit)

	HeatTransferModuleOpt &kernel;

	DiffusionSplit(HeatTransferModuleOpt &kernel): kernel(kernel)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		kernel.gradient.xi.addInputs(kernel.gradient.output, kernel.integration.dND, kernel.material.conductivity, kernel.material.mass);
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
