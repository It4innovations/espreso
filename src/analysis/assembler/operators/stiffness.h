
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_STIFFNESS_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_STIFFNESS_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

#include "config/ecf/material/material.h"
#include "esinfo/meshinfo.h"

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

	Stiffness& operator+=(const size_t rhs)
	{
		dND += rhs; determinant += rhs; conductivity += rhs; xi += rhs; thickness += rhs;
		stiffness += rhs;
		return *this;
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
		ADDMN3M33M3N<nodes>(xi[gpindex] * determinant[gpindex] * weight[gpindex], conductivity.data + 9 * gpindex, dND.data + 3 * nodes * gpindex, stiffness.data);
	}
};

template <class Assembler>
struct HeatStiffness: public ElementOperatorBuilder {
	Assembler &assembler;

	HeatStiffness(Assembler &assembler): ElementOperatorBuilder("HEAT TRANSFER STIFFNESS"), assembler(assembler)
	{

	}

	bool build() override
	{
		if (info::mesh->dimension == 2) {
			assembler.elements.stiffness.addInput(assembler.thickness.gp);
		}
		assembler.elements.stiffness.addInput(assembler.integration.dND);
		assembler.elements.stiffness.addInput(assembler.integration.weight);
		assembler.elements.stiffness.addInput(assembler.integration.jacobiDeterminant);
		assembler.elements.stiffness.addInput(assembler.material.conductivity);
		assembler.elements.stiffness.addInput(assembler.gradient.xi);
		assembler.elements.stiffness.resize();
		assembler.addParameter(assembler.elements.stiffness);
		return true;
	}

	void apply(int interval)
	{
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
		if (mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
			if (info::mesh->dimension == 2) {
				iterate_elements_gps<typename Assembler::NGP>(Stiffness2DHeatIsotropic(assembler.integration.dND, assembler.integration.weight, assembler.integration.jacobiDeterminant, assembler.material.conductivityIsotropic, assembler.gradient.xi, assembler.thickness.gp, assembler.elements.stiffness, interval));
			}
			if (info::mesh->dimension == 3) {
				iterate_elements_gps<typename Assembler::NGP>(Stiffness3DHeatIsotropic(assembler.integration.dND, assembler.integration.weight, assembler.integration.jacobiDeterminant, assembler.material.conductivityIsotropic, assembler.gradient.xi, assembler.thickness.gp, assembler.elements.stiffness, interval));
			}
		} else {
			if (info::mesh->dimension == 2) {
				iterate_elements_gps<typename Assembler::NGP>(Stiffness2DHeat(assembler.integration.dND, assembler.integration.weight, assembler.integration.jacobiDeterminant, assembler.material.conductivity, assembler.gradient.xi, assembler.thickness.gp, assembler.elements.stiffness, interval));
			}
			if (info::mesh->dimension == 3) {
				iterate_elements_gps<typename Assembler::NGP>(Stiffness3DHeat(assembler.integration.dND, assembler.integration.weight, assembler.integration.jacobiDeterminant, assembler.material.conductivity, assembler.gradient.xi, assembler.thickness.gp, assembler.elements.stiffness, interval));
			}
		}
	}
};

//struct Stiffness2DHeatIsotropicSimd: public Stiffness {
//	using Stiffness::Stiffness;
//
//	template<size_t nodes, size_t gps>
//	void operator()(size_t gpindex)
//	{
//		const double * __restrict__ pDND = dND.data;
//		double * __restrict__ pStiffness = stiffness.data;
//
//		SIMD scale =  load(&thickness[gpindex * SIMD::size])
//					* load(&xi[gpindex * SIMD::size])
//					* load(&determinant[gpindex * SIMD::size])
//					* load(&weight[gpindex * SIMD::size])
//					* load(&conductivity[gpindex * SIMD::size]);
//
//		ADDMN2M2NSimd<nodes>(scale, &pDND[2 * nodes * gpindex * SIMD::size], pStiffness);
//	}
//};

//struct HeatStiffnessSimd: public ElementOperatorBuilder {
//	HeatTransferModuleOpt &kernel;
//
//	HeatStiffnessSimd(HeatTransferModuleOpt &kernel): ElementOperatorBuilder("HEAT TRANSFER STIFFNESS"), kernel(kernel)
//	{
//
//	}
//
//	bool build(HeatTransferModuleOpt &kernel) override
//	{
//		if (info::mesh->dimension == 2) {
//			kernel.elementsSimd.stiffness.addInput(kernel.thicknessSimd.gp);
//		}
//		kernel.elementsSimd.stiffness.addInput(kernel.integrationSimd.dND);
//		kernel.elementsSimd.stiffness.addInput(kernel.integrationSimd.weight);
//		kernel.elementsSimd.stiffness.addInput(kernel.integrationSimd.jacobiDeterminant);
//		kernel.elementsSimd.stiffness.addInput(kernel.materialSimd.conductivity);
//		kernel.elementsSimd.stiffness.addInput(kernel.gradientSimd.xi);
//		kernel.elementsSimd.stiffness.resizeAligned(SIMD::size*sizeof(double));
//		kernel.addParameter(kernel.elementsSimd.stiffness);
//		return true;
//	}
//
//	void apply(int interval)
//	{
//		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
//		if (mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
//			if (info::mesh->dimension == 2) {
//				iterate_elements_gps_simd<HeatTransferModuleOpt::NGP>(Stiffness2DHeatIsotropicSimd(kernel.integrationSimd.dND, kernel.integrationSimd.weight, kernel.integrationSimd.jacobiDeterminant, kernel.materialSimd.conductivityIsotropic, kernel.gradientSimd.xi, kernel.thicknessSimd.gp, kernel.elementsSimd.stiffness, interval));
//			}
//		}
//	}
//};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_STIFFNESS_H_ */
