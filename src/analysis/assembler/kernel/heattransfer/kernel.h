
#ifndef SRC_ANALYSIS_ASSEMBLER_KERNEL_HEATTRANSFER_KERNEL_H_
#define SRC_ANALYSIS_ASSEMBLER_KERNEL_HEATTRANSFER_KERNEL_H_

#include "analysis/assembler/module/heattransfer.h"
#include "analysis/assembler/module/assembler.hpp"

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
//	AdvectionKernel<nodes, gps, ndim, etype, Physics> advection(subkernels[interval].advection);
	HeatTransferMatrixKernel<nodes, gps, ndim, model, Physics> K(subkernels.K);
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
//	advection.setActiveness(action);
	K.setActiveness(action);
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
//		if (advection.isactive) {
//			advection.simd(element);
//		}
		if (K.isactive) {
			K.simd(element);
//			if (c == 0) printf("K ");
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

}

#endif /* SRC_ANALYSIS_ASSEMBLER_KERNEL_HEATTRANSFER_KERNEL_H_ */