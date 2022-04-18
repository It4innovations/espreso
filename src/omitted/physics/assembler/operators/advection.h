
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_ADVECTION_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_ADVECTION_H_

#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

struct Advection: public Operator {
	Advection(
			double sigma,
			const ParameterData &tranlationMotion,
			const ParameterData &mass,
			const ParameterData &N,
			const ParameterData &dND,
			const ParameterData &weight,
			const ParameterData &determinant,
			const ParameterData &heatSource,
			const ParameterData &xi,
			ParameterData &conductivity,
			ParameterData &stiffness,
			ParameterData &rhs,
			int interval)
	: Operator(interval,
			conductivity.isconst[interval] && stiffness.isconst[interval] && rhs.isconst[interval],
			conductivity.update[interval] || stiffness.update[interval] || rhs.update[interval]),
	  sigma(sigma),
	  translationMotion(tranlationMotion, interval),
	  mass(mass, interval),
	  N(N, interval),
	  dND(dND, interval),
	  weight(weight, interval, 0),
	  determinant(determinant, interval),
	  heatSource(heatSource, interval),
	  xi(xi, interval),
	  conductivity(conductivity, interval),
	  stiffness(stiffness, interval),
	  rhs(rhs, interval)
	{

	}

	double sigma;
	InputParameterIterator translationMotion, mass, N, dND, weight, determinant, heatSource, xi;
	OutputParameterIterator conductivity, stiffness, rhs;

	void operator++()
	{
		++translationMotion; ++mass; ++N; ++dND; ++determinant, ++xi;
		++conductivity; ++stiffness; ++rhs;
	}
};

template<bool Isotropic, bool HeatSource, bool CAU, bool DiffusionSplit>
struct Advection2D: public Advection {
	using Advection::Advection;

	template<int size>
	double norm(double *v)
	{
		double norm = 0;
		for (int n = 0; n < size; ++n) {
			norm += v[n] * v[n];
		}
		return std::sqrt(norm);
	}

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		double u[2] = { translationMotion[2 * gpindex] * mass[gpindex], translationMotion[2 * gpindex + 1] * mass[gpindex] };
		double scale = 0;
		if (u[0] || u[1]) {
			double b[nodes];
			M12M2N<nodes>(u, dND.data + 2 * gpindex * nodes, b);

			double unorm = std::sqrt(u[0] * u[0] + u[1] * u[1]);
			double bnorm = norm<nodes>(b);
			double h = 2 * unorm / bnorm;
			double P = h * unorm / (2 * conductivity[0]);
			double T = std::max(.0, 1 - 1 / P);
			scale = h * T / (2 * unorm);

			ADDMN1M1N<nodes>(xi.data[gpindex] * weight[gpindex] * determinant[gpindex], N.data, b, stiffness.data);
			ADDMN1M1N<nodes>(scale * xi.data[gpindex] * weight[gpindex] * determinant[gpindex], b, stiffness.data);

			if (Isotropic) {
				conductivity[0] += sigma * h * unorm;
			} else {
				conductivity[0] += sigma * h * unorm;
				conductivity[3] += sigma * h * unorm;
			}

			if (HeatSource) {
				for (int n = 0; n < nodes; ++n) {
					rhs.data[n] += determinant[gpindex] * weight[gpindex] * h * T * b[n] * heatSource[n] / (2 * unorm);
				}
			}

			if (CAU) {
				double dNDNorm = norm<2 * nodes>(dND);
				double v[2] = { 0, 0 };
				double r[nodes];
				if (dNDNorm >= 1e-12) {
					for (int n = 0; n < nodes; n++) {
						r[n] = b[n] - heatSource[n];
					}
					M1NMN2(1 / (dNDNorm * dNDNorm), r, dND.data + 2 * gpindex * nodes, v);
				}

				double c[nodes];
				M12M2N<nodes>(v, dND.data + 2 * gpindex * nodes, c);
				double vnorm = std::sqrt(v[0] * v[0] + v[1] * v[1]);
				double cnorm = norm<nodes>(c);
				double hc = 2 * vnorm / cnorm;
				double cenorm = 0;
				if (Isotropic) {
					cenorm = conductivity[0];
				} else {
					cenorm = norm<4>(conductivity.data);
				}
				double Pc = hc * vnorm / (2 * cenorm);
				double Tc = std::max(.0, 1 - 1 / Pc);

				double c1 = dNDNorm >= 1e-12 ? c1 = norm<nodes>(r) : bnorm;
				double c2 = T * h != 0 ? Tc * hc / (T * h) : 0;
				double C = c2 / unorm < c2 ? T * h * c1 * (c2 - c1 / unorm) / 2 : 0;
				if (c1 / unorm < c2) {
					double C = T * h * c1 * (c2 - c1 / unorm) / 2;
					ADDMN2M2N<nodes>(C * xi.data[gpindex] * weight[gpindex] * determinant[gpindex], dND.data + 2 * gpindex * nodes, stiffness.data);
				}
			}
		}
	}
};

struct Advection3DIsotropic: public Advection {
	using Advection::Advection;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{

	}
};

struct TranslationMotion: public ElementOperatorBuilder {
	HeatTransferModuleOpt &kernel;

	TranslationMotion(HeatTransferModuleOpt &kernel): ElementOperatorBuilder("TRANSLATION MOTION"), kernel(kernel)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		kernel.material.conductivityIsotropic.addInput(kernel.translationMotions.gp);
		kernel.material.conductivityIsotropic.addInput(kernel.material.mass);
		kernel.material.conductivityIsotropic.resize();

		kernel.material.conductivity.addInput(kernel.translationMotions.gp);
		kernel.material.conductivity.addInput(kernel.material.mass);
		kernel.material.conductivity.resize();

		kernel.translationMotions.stiffness.addInput(kernel.material.mass);
		kernel.translationMotions.stiffness.addInput(kernel.translationMotions.gp);
		kernel.translationMotions.stiffness.resize();
		kernel.addParameter(kernel.translationMotions.stiffness);
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
