
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_

#include "copy.h"
#include "coordinatesystem.h"
#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"
#include "physics/assembler/math.hpp"

#include "esinfo/ecfinfo.h"

namespace espreso {

struct CopyConductivity: public Operator {
	GET_NAME(CopyConductivity)
	CopyConductivity(
			const ParameterData &input,
			ParameterData &output,
			int interval)
	: Operator(interval, output.isconst[interval], Link(interval).inputs(input).outputs(output)),
	  input(input, interval, input.increment(interval)),
	  output(output, interval, ndim * ndim * egps)
	{ }

	InputParameterIterator input;
	OutputParameterIterator output;

	void operator++()
	{
		++input;
		++output;
	}
};

struct CopyDiagonal2DConductivity: public CopyConductivity {
	GET_NAME(CopyDiagonal2DConductivity)
	using CopyConductivity::CopyConductivity;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		int igp = 2 * gpindex;
		int ogp = 4 * gpindex;
		output.data[ogp + 0] = input.data[igp + 0]; output.data[ogp + 1] = 0;
		output.data[ogp + 2] = 0;                   output.data[ogp + 3] = input.data[igp + 1];
	}
};

struct CopyDiagonal3DConductivity: public CopyConductivity {
	GET_NAME(CopyDiagonal3DConductivity)
	using CopyConductivity::CopyConductivity;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		int igp = 3 * gpindex;
		int ogp = 9 * gpindex;
		output.data[ogp + 0] = input.data[igp + 0]; output.data[ogp + 1] = 0;                   output.data[ogp + 2] = 0;
		output.data[ogp + 3] = 0;                   output.data[ogp + 4] = input.data[igp + 1]; output.data[ogp + 5] = 0;
		output.data[ogp + 6] = 0;                   output.data[ogp + 7] = 0;                   output.data[ogp + 8] = input.data[igp + 2];
	}
};

struct CopySymmetric2DConductivity: public CopyConductivity {
	GET_NAME(CopySymmetric2DConductivity)
	using CopyConductivity::CopyConductivity;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		int igp = 3 * gpindex;
		int ogp = 4 * gpindex;
		output.data[ogp + 0] = input.data[igp + 0]; output.data[ogp + 1] = input.data[igp + 2];
		output.data[ogp + 2] = input.data[igp + 2]; output.data[ogp + 3] = input.data[igp + 1];
	}
};

struct CopySymmetric3DConductivity: public CopyConductivity {
	GET_NAME(CopySymmetric3DConductivity)
	using CopyConductivity::CopyConductivity;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		int igp = 6 * gpindex;
		int ogp = 9 * gpindex;
		output.data[ogp + 0] = input.data[igp + 0]; output.data[ogp + 1] = input.data[igp + 3]; output.data[ogp + 2] = input.data[igp + 5];
		output.data[ogp + 3] = 0;                   output.data[ogp + 4] = input.data[igp + 1]; output.data[ogp + 5] = input.data[igp + 4];
		output.data[ogp + 6] = 0;                   output.data[ogp + 7] = 0;                   output.data[ogp + 8] = input.data[igp + 2];
	}
};

struct CopyAnisotropic2DConductivity: public CopyConductivity {
	GET_NAME(CopyAnisotropic2DConductivity)
	using CopyConductivity::CopyConductivity;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		int igp = 4 * gpindex;
		int ogp = 4 * gpindex;
		output.data[ogp + 0] = input.data[igp + 0]; output.data[ogp + 1] = input.data[igp + 2];
		output.data[ogp + 2] = input.data[igp + 3]; output.data[ogp + 3] = input.data[igp + 1];
	}
};


struct CopyAnisotropic3DConductivity: public CopyConductivity {
	GET_NAME(CopyAnisotropic3DConductivity)
	using CopyConductivity::CopyConductivity;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		int igp = 9 * gpindex;
		int ogp = 9 * gpindex;
		output.data[ogp + 0] = input.data[igp + 0]; output.data[ogp + 1] = input.data[igp + 3]; output.data[ogp + 2] = input.data[igp + 5];
		output.data[ogp + 3] = input.data[igp + 6]; output.data[ogp + 4] = input.data[igp + 1]; output.data[ogp + 5] = input.data[igp + 4];
		output.data[ogp + 6] = input.data[igp + 8]; output.data[ogp + 7] = input.data[igp + 7]; output.data[ogp + 8] = input.data[igp + 2];
	}
};

struct ThermalConductivity: public ElementOperatorBuilder {
	GET_NAME(ThermalConductivity)

	HeatTransferModuleOpt &kernel;

	ThermalConductivity(HeatTransferModuleOpt &kernel): kernel(kernel)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[i].material];
			switch (mat->coordinate_system.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN: break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL: kernel.material.conductivity.addInput(i, info::mesh->nodes->coordinates); break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: kernel.material.conductivity.addInput(i, info::mesh->nodes->coordinates); break;
				break;
			}
		}

		kernel.material.conductivityIsotropic.addInputs(kernel.material.model.isotropic);
		if (info::mesh->dimension == 2) {
			kernel.material.conductivity.addInputs(kernel.material.model.diagonal, kernel.material.model.symmetric2D, kernel.material.model.anisotropic, kernel.cooSystem.cartesian2D, kernel.cooSystem.spherical, kernel.cooSystem.cylindric);
		}
		if (info::mesh->dimension == 3) {
			kernel.material.conductivity.addInputs(kernel.material.model.diagonal, kernel.material.model.symmetric3D, kernel.material.model.anisotropic, kernel.cooSystem.cartesian3D, kernel.cooSystem.spherical, kernel.cooSystem.cylindric);
		}
		return true;
	}

	void apply(int interval)
	{
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
		if (info::mesh->dimension == 2) {
			switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: CopyElementParameters(kernel.material.model.isotropic, kernel.material.conductivityIsotropic).apply(interval); break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL: iterate_elements_gps<HeatTransferModuleOpt>(CopyDiagonal2DConductivity(kernel.material.model.diagonal, kernel.material.conductivity, interval)); break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC: iterate_elements_gps<HeatTransferModuleOpt>(CopySymmetric2DConductivity(kernel.material.model.symmetric2D, kernel.material.conductivity, interval)); break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: iterate_elements_gps<HeatTransferModuleOpt>(CopyAnisotropic2DConductivity(kernel.material.model.anisotropic, kernel.material.conductivity, interval)); break;
				break;
			}

			switch (mat->coordinate_system.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN: iterate_elements_gps<HeatTransferModuleOpt>(Cartesian2DCoordinateSystem(kernel.coords.gp, kernel.cooSystem.cartesian2D, kernel.material.conductivity, interval)); break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL: break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: iterate_elements_gps<HeatTransferModuleOpt>(Cylindrical2DCoordinateSystem(kernel.coords.gp, kernel.cooSystem.cylindric, kernel.material.conductivity, interval)); break;
				break;
			}
		}
		if (info::mesh->dimension == 3) {
			switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: CopyElementParameters(kernel.material.model.isotropic, kernel.material.conductivityIsotropic).apply(interval); break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL: iterate_elements_gps<HeatTransferModuleOpt>(CopyDiagonal3DConductivity(kernel.material.model.diagonal, kernel.material.conductivity, interval)); break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC: iterate_elements_gps<HeatTransferModuleOpt>(CopyDiagonal3DConductivity(kernel.material.model.symmetric3D, kernel.material.conductivity, interval)); break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: iterate_elements_gps<HeatTransferModuleOpt>(CopyDiagonal3DConductivity(kernel.material.model.anisotropic, kernel.material.conductivity, interval)); break;
				break;
			}

			switch (mat->coordinate_system.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN: iterate_elements_gps<HeatTransferModuleOpt>(Cartesian3DCoordinateSystem(kernel.coords.gp, kernel.cooSystem.cartesian3D, kernel.material.conductivity, interval)); break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL: iterate_elements_gps<HeatTransferModuleOpt>(Spherical3DCoordinateSystem(kernel.coords.gp, kernel.cooSystem.spherical, kernel.material.conductivity, interval)); break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: iterate_elements_gps<HeatTransferModuleOpt>(Cylindrical3DCoordinateSystem(kernel.coords.gp, kernel.cooSystem.cylindric, kernel.material.conductivity, interval)); break;
				break;
			}
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_ */
