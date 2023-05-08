
#ifndef SRC_ANALYSIS_ASSEMBLER_KERNEL_HEATTRANSFER_KERNEL_H_
#define SRC_ANALYSIS_ASSEMBLER_KERNEL_HEATTRANSFER_KERNEL_H_

#include "analysis/assembler/module/heattransfer.h"
#include "analysis/assembler/module/assembler.hpp"
#include "analysis/assembler/subkernel/heattransfer/heatsource.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, enum ThermalConductivityConfiguration::MODEL ecfmodel, enum ThermalConductivityConfiguration::MODEL model>
void compute(const HeatTransfer::SubKernels &subkernels, Assembler::Action action)
{
	typedef HeatTransferElementDescriptor<nodes, gps, ndim, edim, ecfmodel, model> Physics;
	typename Physics::Element element;

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, gps, ndim, Physics> coordinates(subkernels.coordinates);
	ThicknessToNodes<nodes, ndim, Physics> thickness(subkernels.thickness);
	TemperatureKernel<nodes, gps, Physics> temperature(subkernels.temperature);
	IntegrationKernel<nodes, gps, ndim, edim, Physics> integration(subkernels.integration);
	HeatTransferCoordinateSystemKernel<gps, ndim, ecfmodel, model, Physics> coosystem(subkernels.coosystem);
	HeatSourceKernel<nodes, gps, ndim, Physics> heatSource(subkernels.heatSource);
	AdvectionKernel<nodes, gps, ndim, model, Physics> advection(subkernels.advection);
	HeatTransferMatrixKernel<nodes, gps, ndim, model, Physics> K(subkernels.K), M(subkernels.M);
	AdvectionMatrix<nodes, gps, ndim, model, Physics> advK(subkernels.advection);
	TemperatureGradientKernel<nodes, gps, ndim, Physics> gradient(subkernels.gradient);
	TemperatureFluxKernel<nodes, gps, ndim, model, Physics> flux(subkernels.flux);

	std::vector<ExternalGPsExpression<gps, Physics>*> nonconst;
	for (size_t i = 0; i < subkernels.expressions.size(); ++i) {
		if (subkernels.expressions[i]->evaluator->isConst()) {
			dynamic_cast<ExternalGPsExpression<gps, Physics>*>(subkernels.expressions[i])->simd(element);
		} else {
			nonconst.push_back(dynamic_cast<ExternalGPsExpression<gps, Physics>*>(subkernels.expressions[i]));
		}
	}

	basis.simd(element);
	if (coosystem.isactive) {
		coosystem.simd(element);
	}

	thickness.setActiveness(action);
	temperature.setActiveness(action);
	coosystem.setActiveness(action);
	heatSource.setActiveness(action);
	advection.setActiveness(action);
	K.setActiveness(action);
	M.setActiveness(action);
	advK.setActiveness(action);
	gradient.setActiveness(action);
	flux.setActiveness(action);

//	printf("sub-kernels: ");
	for (esint c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
//		if (c == 0) printf("coordinates ");
		if (temperature.isactive) {
			temperature.simd(element);
//			if (c == 0) printf("temp ");
		}
		if (thickness.isactive) {
			thickness.simd(element);
//			if (c == 0) printf("thickness ");
		}
		integration.simd(element);
//		if (c == 0) printf("integrate ");
		if (coosystem.isactive) {
			coosystem.simd(element);
//			if (c == 0) printf("coosystem ");
		}
		if (heatSource.isactive) {
			heatSource.simd(element);
//			if (c == 0) printf("heatSource ");
		}
		if (advection.isactive) {
			advection.simd(element);
//			if (c == 0) printf("advection ");
		}
		if (K.isactive) {
			K.simd(element);
//			if (c == 0) printf("K ");
		}
		if (advK.isactive) {
			advK.simd(element);
//			if (c == 0) printf("advK ");
		}
		if (gradient.isactive) {
			gradient.simd(element);
//			if (c == 0) printf("gradient ");
		}
		if (flux.isactive) {
			flux.simd(element);
//			if (c == 0) printf("flux ");
		}
	}
//	printf("\n");
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runConductivity(HeatTransfer::SubKernels &subkernels, Assembler::Action action)
{
	switch (subkernels.conductivity.conductivity->model) {
	case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
		compute<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::ISOTROPIC, ThermalConductivityConfiguration::MODEL::ISOTROPIC>(subkernels, action);
		break;
	case ThermalConductivityConfiguration::MODEL::DIAGONAL:
		if (subkernels.coosystem.rotated) {
			compute<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::DIAGONAL, ThermalConductivityConfiguration::MODEL::SYMMETRIC>(subkernels, action);
		} else {
			compute<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::DIAGONAL, ThermalConductivityConfiguration::MODEL::DIAGONAL>(subkernels, action);
		}
		break;
	case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
		compute<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::SYMMETRIC, ThermalConductivityConfiguration::MODEL::SYMMETRIC>(subkernels, action);
		break;
	case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
		compute<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::ANISOTROPIC, ThermalConductivityConfiguration::MODEL::ANISOTROPIC>(subkernels, action);
		break;
	}
}

}

#endif /* SRC_ANALYSIS_ASSEMBLER_KERNEL_HEATTRANSFER_KERNEL_H_ */
