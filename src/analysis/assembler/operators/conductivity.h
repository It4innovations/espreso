
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_

#include "copy.h"
#include "coordinatesystem.h"
#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

#include "esinfo/ecfinfo.h"

namespace espreso {

template <class Assembler>
struct ThermalConductivity: public ElementOperatorBuilder {

	struct CopyConductivity: public Operator {
		CopyConductivity(
				const ParameterData &input,
				ParameterData &output,
				int interval)
		: Operator(interval, output.isconst[interval], output.update[interval]),
		  input(input, interval),
		  output(output, interval)
		{

		}

		InputParameterIterator input;
		OutputParameterIterator output;

		void operator++()
		{
			++input;
			++output;
		}
	};

	struct CopyDiagonal2DConductivity: public CopyConductivity {
		using CopyConductivity::CopyConductivity;

		template<int nodes, int gps>
		void operator()(int gpindex)
		{
			int igp = 2 * gpindex;
			int ogp = 4 * gpindex;
			this->output.data[ogp + 0] = this->input.data[igp + 0]; this->output.data[ogp + 1] = 0;
			this->output.data[ogp + 2] = 0;                         this->output.data[ogp + 3] = this->input.data[igp + 1];
		}
	};

	struct CopyDiagonal3DConductivity: public CopyConductivity {
		using CopyConductivity::CopyConductivity;

		template<int nodes, int gps>
		void operator()(int gpindex)
		{
			int igp = 3 * gpindex;
			int ogp = 9 * gpindex;
			this->output.data[ogp + 0] = this->input.data[igp + 0]; this->output.data[ogp + 1] = 0;                         this->output.data[ogp + 2] = 0;
			this->output.data[ogp + 3] = 0;                         this->output.data[ogp + 4] = this->input.data[igp + 1]; this->output.data[ogp + 5] = 0;
			this->output.data[ogp + 6] = 0;                         this->output.data[ogp + 7] = 0;                         this->output.data[ogp + 8] = this->input.data[igp + 2];
		}
	};

	struct CopySymmetric2DConductivity: public CopyConductivity {
		using CopyConductivity::CopyConductivity;

		template<int nodes, int gps>
		void operator()(int gpindex)
		{
			int igp = 3 * gpindex;
			int ogp = 4 * gpindex;
			this->output.data[ogp + 0] = this->input.data[igp + 0]; this->output.data[ogp + 1] = this->input.data[igp + 2];
			this->output.data[ogp + 2] = this->input.data[igp + 2]; this->output.data[ogp + 3] = this->input.data[igp + 1];
		}
	};

	struct CopySymmetric3DConductivity: public CopyConductivity {
		using CopyConductivity::CopyConductivity;

		template<int nodes, int gps>
		void operator()(int gpindex)
		{
			int igp = 6 * gpindex;
			int ogp = 9 * gpindex;
			this->output.data[ogp + 0] = this->input.data[igp + 0]; this->output.data[ogp + 1] = this->input.data[igp + 3]; this->output.data[ogp + 2] = this->input.data[igp + 5];
			this->output.data[ogp + 3] = this->input.data[igp + 3]; this->output.data[ogp + 4] = this->input.data[igp + 1]; this->output.data[ogp + 5] = this->input.data[igp + 4];
			this->output.data[ogp + 6] = this->input.data[igp + 5]; this->output.data[ogp + 7] = this->input.data[igp + 4]; this->output.data[ogp + 8] = this->input.data[igp + 2];
		}
	};

	struct CopyAnisotropic2DConductivity: public CopyConductivity {
		using CopyConductivity::CopyConductivity;

		template<int nodes, int gps>
		void operator()(int gpindex)
		{
			int igp = 4 * gpindex;
			int ogp = 4 * gpindex;
			this->output.data[ogp + 0] = this->input.data[igp + 0]; this->output.data[ogp + 1] = this->input.data[igp + 2];
			this->output.data[ogp + 2] = this->input.data[igp + 3]; this->output.data[ogp + 3] = this->input.data[igp + 1];
		}
	};


	struct CopyAnisotropic3DConductivity: public CopyConductivity {
		using CopyConductivity::CopyConductivity;

		template<int nodes, int gps>
		void operator()(int gpindex)
		{
			int igp = 9 * gpindex;
			int ogp = 9 * gpindex;
			this->output.data[ogp + 0] = this->input.data[igp + 0]; this->output.data[ogp + 1] = this->input.data[igp + 3]; this->output.data[ogp + 2] = this->input.data[igp + 5];
			this->output.data[ogp + 3] = this->input.data[igp + 6]; this->output.data[ogp + 4] = this->input.data[igp + 1]; this->output.data[ogp + 5] = this->input.data[igp + 4];
			this->output.data[ogp + 6] = this->input.data[igp + 8]; this->output.data[ogp + 7] = this->input.data[igp + 7]; this->output.data[ogp + 8] = this->input.data[igp + 2];
		}
	};

	Assembler &assembler;

	ThermalConductivity(Assembler &assembler): ElementOperatorBuilder("THERMAL CONDUCTIVITY MATRIX"), assembler(assembler)
	{

	}

	bool build() override
	{
		for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[i].material];
			switch (mat->coordinate_system.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN: break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL: assembler.material.conductivity.addInput(i, info::mesh->nodes->coordinates); break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: assembler.material.conductivity.addInput(i, info::mesh->nodes->coordinates); break;
				break;
			}
		}

		assembler.material.conductivityIsotropic.addInput(assembler.material.model.isotropic);
		assembler.material.conductivity.addInput(assembler.material.model.diagonal);
		assembler.material.conductivity.addInput(assembler.material.model.anisotropic);
		assembler.material.conductivity.addInput(assembler.cooSystem.spherical);
		assembler.material.conductivity.addInput(assembler.cooSystem.cylindric);
		if (info::mesh->dimension == 2) {
			assembler.material.conductivity.addInput(assembler.material.model.symmetric2D);
			assembler.material.conductivity.addInput(assembler.cooSystem.cartesian2D);
		}
		if (info::mesh->dimension == 3) {
			assembler.material.conductivity.addInput(assembler.material.model.symmetric3D);
			assembler.material.conductivity.addInput(assembler.cooSystem.cartesian3D);
		}

		assembler.material.conductivityIsotropic.resize();
		assembler.material.conductivity.resize();

		assembler.addParameter(assembler.material.conductivityIsotropic);
		assembler.addParameter(assembler.material.conductivity);
		return true;
	}

	void apply(int interval)
	{
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
		if (info::mesh->dimension == 2) {
			switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: CopyElementParameters<Assembler>(assembler, assembler.material.model.isotropic, assembler.material.conductivityIsotropic, "COPY ISOTROPIC CONDUCTIVITY").apply(interval); return;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL: iterate_elements_gps<typename Assembler::NGP>(CopyDiagonal2DConductivity(assembler.material.model.diagonal, assembler.material.conductivity, interval)); break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC: iterate_elements_gps<typename Assembler::NGP>(CopySymmetric2DConductivity(assembler.material.model.symmetric2D, assembler.material.conductivity, interval)); break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: iterate_elements_gps<typename Assembler::NGP>(CopyAnisotropic2DConductivity(assembler.material.model.anisotropic, assembler.material.conductivity, interval)); break;
				break;
			}

			switch (mat->coordinate_system.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN: iterate_elements_gps<typename Assembler::NGP>(Cartesian2DCoordinateSystem(assembler.coords.gp, assembler.cooSystem.cartesian2D, assembler.material.conductivity, interval)); break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL: break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: iterate_elements_gps<typename Assembler::NGP>(Cylindrical2DCoordinateSystem(assembler.coords.gp, assembler.cooSystem.cylindric, assembler.material.conductivity, interval)); break;
				break;
			}
		}
		if (info::mesh->dimension == 3) {
			switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: CopyElementParameters<Assembler>(assembler, assembler.material.model.isotropic, assembler.material.conductivityIsotropic, "COPY ISOTROPIC CONDUCTIVITY").apply(interval); return;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL: iterate_elements_gps<typename Assembler::NGP>(CopyDiagonal3DConductivity(assembler.material.model.diagonal, assembler.material.conductivity, interval)); break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC: iterate_elements_gps<typename Assembler::NGP>(CopySymmetric3DConductivity(assembler.material.model.symmetric3D, assembler.material.conductivity, interval)); break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: iterate_elements_gps<typename Assembler::NGP>(CopyAnisotropic3DConductivity(assembler.material.model.anisotropic, assembler.material.conductivity, interval)); break;
				break;
			}

			switch (mat->coordinate_system.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN: iterate_elements_gps<typename Assembler::NGP>(Cartesian3DCoordinateSystem(assembler.coords.gp, assembler.cooSystem.cartesian3D, assembler.material.conductivity, interval)); break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL: iterate_elements_gps<typename Assembler::NGP>(Spherical3DCoordinateSystem(assembler.coords.gp, assembler.cooSystem.spherical, assembler.material.conductivity, interval)); break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: iterate_elements_gps<typename Assembler::NGP>(Cylindrical3DCoordinateSystem(assembler.coords.gp, assembler.cooSystem.cylindric, assembler.material.conductivity, interval)); break;
				break;
			}
		}
	}
};

//struct ThermalConductivitySimd: public ElementOperatorBuilder {
//
//	struct CopyConductivity: public Operator {
//		CopyConductivity(
//				const ParameterData &input,
//				ParameterData &output,
//				int interval)
//		: Operator(interval, output.isconst[interval], output.update[interval]),
//		  input(input, interval),
//		  output(output, interval)
//		{
//
//		}
//
//		InputParameterIterator input;
//		OutputParameterIterator output;
//
//		void operator++()
//		{
//			++input;
//			++output;
//		}
//	};
//
//	HeatTransferModuleOpt &kernel;
//
//	ThermalConductivitySimd(HeatTransferModuleOpt &kernel): ElementOperatorBuilder("THERMAL CONDUCTIVITY MATRIX"), kernel(kernel)
//	{
//
//	}
//
//	bool build(HeatTransferModuleOpt &kernel) override
//	{
//		for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
//			const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[i].material];
//			switch (mat->coordinate_system.type) {
//			case CoordinateSystemConfiguration::TYPE::CARTESIAN: break;
//				break;
//			}
//		}
//
//		kernel.materialSimd.conductivityIsotropic.addInput(kernel.material.model.isotropic);
//		kernel.materialSimd.conductivityIsotropic.resizeAligned(SIMD::size*sizeof(double));
//		kernel.addParameter(kernel.materialSimd.conductivityIsotropic);
//		return true;
//	}
//
//	void apply(int interval)
//	{
//		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
//		if (info::mesh->dimension == 2) {
//			switch (mat->thermal_conductivity.model) {
//			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: CopyElementParameters(kernel.materialSimd.model.isotropic, kernel.materialSimd.conductivityIsotropic, "COPY ISOTROPIC CONDUCTIVITY").apply(interval); return;
//			}
//		}
//		if (info::mesh->dimension == 3) {
//			switch (mat->thermal_conductivity.model) {
//			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: CopyElementParameters(kernel.materialSimd.model.isotropic, kernel.materialSimd.conductivityIsotropic, "COPY ISOTROPIC CONDUCTIVITY").apply(interval); return;
//			}
//		}
//	}
//};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_ */
