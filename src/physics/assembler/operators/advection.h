
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_ADVECTION_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_ADVECTION_H_

#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

struct Advection: public Operator {
	Advection(
			double sigma,
			ParameterData &tranlationMotion,
			ParameterData &mass,
			ParameterData &dND,
			const ParameterData &weight,
			const ParameterData &determinant,
			ParameterData &conductivity,
			ParameterData &stiffness,
			int interval)
	: Operator(interval, tranlationMotion.isconst[interval], Link(interval).inputs(tranlationMotion, mass, dND).outputs(stiffness).self(conductivity)),
	  sigma(sigma),
	  translationMotion(tranlationMotion, interval, egps),
	  mass(mass, interval, egps),
	  dND(dND, interval, dND.size),
	  weight(weight, interval, 0),
	  determinant(determinant, interval, egps),
	  conductivity(conductivity, interval, conductivity.size),
	  stiffness(stiffness, interval, enodes * enodes)
	{

	}

	double sigma;
	InputParameterIterator translationMotion, mass, dND, weight, determinant;
	OutputParameterIterator conductivity, stiffness;

	void operator++()
	{
		++translationMotion; ++mass, ++dND; ++determinant;
		++conductivity; ++stiffness;
	}
};

struct Advection2DIsotropic: public Advection {
	GET_NAME(Advection2DIsotropic)
	using Advection::Advection;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		double u[2] = { translationMotion[0] * mass[0], translationMotion[1] * mass[1] };
		double konst = 0;

		double b[nodes];
		M12M2N<nodes>(u, dND.data + 2 * gpindex * nodes, b);

		if (u[0] || u[1]) {
			double unorm = std::sqrt(u[0] * u[0] + u[1] * u[1]);
			double bnorm = 0;
			for (int n = 0; n < nodes; ++n) {
				bnorm += b[n] * b[n];
			}
			bnorm = std::sqrt(bnorm);
			double h = 2 * unorm / bnorm;
			double P = h * unorm / (2 * conductivity[0]);
			double T = std::max(.0, 1 - 1 / P);
			konst = h * T / (2 * unorm);
			conductivity[0] += sigma * h * unorm;

			ADDMN2M2N<nodes>(determinant[gpindex] * weight[gpindex] * conductivity[gpindex], dND.data + 2 * nodes * gpindex, stiffness.data);
		}
	}
};

struct Advection3DIsotropic: public Advection {
	GET_NAME(Advection3DIsotropic)
	using Advection::Advection;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{

	}
};

struct TranslationMotion: public ElementOperatorBuilder {
	GET_NAME(Gradient)

	HeatTransferModuleOpt &kernel;

	TranslationMotion(HeatTransferModuleOpt &kernel): kernel(kernel)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		kernel.translationMotions.stiffness.addInputs(kernel.material.mass, kernel.translationMotions.gp);
		return true;
	}

	void apply(int interval)
	{
//		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
//		if (mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
//			if (info::mesh->dimension == 2) {
//				iterate_elements_gps<HeatTransferModuleOpt>(Advection2DIsotropic(kernel.translationMotions.sigma, kernel.translationMotions.gp, kernel.material.mass, kernel.integration.dND, kernel.material.conductivityIsotropic, kernel.translationMotions.stiffness, interval));
//			}
//			if (info::mesh->dimension == 3) {
//				iterate_elements_gps<HeatTransferModuleOpt>(Advection3DIsotropic(kernel.translationMotions.sigma, kernel.translationMotions.gp, kernel.material.mass, kernel.integration.dND, kernel.material.conductivityIsotropic, kernel.translationMotions.stiffness, interval));
//			}
//		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_ADVECTION_H_ */
