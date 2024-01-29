
#ifndef SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_KERNEL_H_
#define SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_KERNEL_H_

#include "element.h"
#include "operators.h"
#include "mesh/element.h"

#include <iostream>

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void setElementKernel(HeatTransferElementOperators &subkernels, SubKernel::Action action)
{
	typedef HeatTransferElement<nodes, gps, ndim, edim> Element; Element element;

	if constexpr(ndim == 2) {
		if (subkernels.thickness.expression) {
			subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
					subkernels.thickness.expression->evaluator,
					[] (Element &element, size_t &n, size_t &s, double value) { element.thickness.node[n][s] = value; }));
		}

		switch (subkernels.conductivity.coordinateSystem->type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN:
			if (subkernels.conductivity.coordinateSystem->rotation.z.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.conductivity.coordinateSystem->rotation.z.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
			}
			break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			if (subkernels.conductivity.coordinateSystem->center.x.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.conductivity.coordinateSystem->center.x.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
			}
			if (subkernels.conductivity.coordinateSystem->center.y.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.conductivity.coordinateSystem->center.y.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
			}
			break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL:
			break;
		}

		Evaluator *kxx = subkernels.conductivity.conductivity->values.get(0, 0).evaluator;
		Evaluator *kxy = subkernels.conductivity.conductivity->values.get(0, 1).evaluator;
		Evaluator *kyx = subkernels.conductivity.conductivity->values.get(1, 0).evaluator;
		Evaluator *kyy = subkernels.conductivity.conductivity->values.get(1, 1).evaluator;
		if (subkernels.conductivity.conductivity->model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				kxx, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[0][s] = element.ecf.conductivity[3][s] = value; }));
		} else {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				kxx, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[0][s] = value; }));
		}
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
			kxy, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[1][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
			kyx, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[2][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
			kyy, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[3][s] = value; }));

		if (subkernels.advection.expressionVector) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.advection.expressionVector->x.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[0][s] = value; }));
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.advection.expressionVector->y.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[1][s] = value; }));
		}
	}

	if constexpr(ndim == 3) {
		switch (subkernels.conductivity.coordinateSystem->type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN:
			if (subkernels.conductivity.coordinateSystem->rotation.x.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.conductivity.coordinateSystem->rotation.x.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
			}
			if (subkernels.conductivity.coordinateSystem->rotation.y.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.conductivity.coordinateSystem->rotation.y.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
			}
			if (subkernels.conductivity.coordinateSystem->rotation.z.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.conductivity.coordinateSystem->rotation.z.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[2][s] = value; }));
			}
			break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			if (subkernels.conductivity.coordinateSystem->center.x.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.conductivity.coordinateSystem->center.x.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
			}
			if (subkernels.conductivity.coordinateSystem->center.y.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.conductivity.coordinateSystem->center.y.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
			}
			break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL:
			if (subkernels.conductivity.coordinateSystem->center.x.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.conductivity.coordinateSystem->center.x.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
			}
			if (subkernels.conductivity.coordinateSystem->center.y.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.conductivity.coordinateSystem->center.y.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
			}
			if (subkernels.conductivity.coordinateSystem->center.z.isset) {
				subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.conductivity.coordinateSystem->center.z.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[2][s] = value; }));
			}
			break;
		}

		Evaluator *kxx = subkernels.conductivity.conductivity->values.get(0, 0).evaluator;
		Evaluator *kxy = subkernels.conductivity.conductivity->values.get(0, 1).evaluator;
		Evaluator *kxz = subkernels.conductivity.conductivity->values.get(0, 2).evaluator;
		Evaluator *kyx = subkernels.conductivity.conductivity->values.get(1, 0).evaluator;
		Evaluator *kyy = subkernels.conductivity.conductivity->values.get(1, 1).evaluator;
		Evaluator *kyz = subkernels.conductivity.conductivity->values.get(1, 2).evaluator;
		Evaluator *kzx = subkernels.conductivity.conductivity->values.get(2, 0).evaluator;
		Evaluator *kzy = subkernels.conductivity.conductivity->values.get(2, 1).evaluator;
		Evaluator *kzz = subkernels.conductivity.conductivity->values.get(2, 2).evaluator;

		if (subkernels.conductivity.conductivity->model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				kxx, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[0][s] = element.ecf.conductivity[4][s] = element.ecf.conductivity[9][s] = value; }));
		} else {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				kxx, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[0][s] = value; }));
		}
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
			kxy, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[1][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
			kxz, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[2][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
			kyx, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[3][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
			kyy, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[4][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
			kyz, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[5][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
			kzx, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[6][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
			kzy, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[7][s] = value; }));
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
			kzz, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[8][s] = value; }));

		if (subkernels.advection.expressionVector) {
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.advection.expressionVector->x.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[0][s] = value; }));
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.advection.expressionVector->y.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[1][s] = value; }));
			subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
					subkernels.advection.expressionVector->z.evaluator,
					[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[2][s] = value; }));
		}
	}

	if (subkernels.initTemperature.expression) {
		subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
				subkernels.initTemperature.expression->evaluator,
				[] (Element &element, size_t &n, size_t &s, double value) { element.temperature.initial[n][s] = value; }));
	}

	if (subkernels.heatSource.expression) {
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.heatSource.expression->evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.heatSource[s] = value; }));
	}

	if (subkernels.material.configuration->density.isset) {
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.material.configuration->density.evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.density[s] = value; }));
	}
	if (subkernels.material.configuration->heat_capacity.isset) {
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.material.configuration->heat_capacity.evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.heatCapacity[s] = value; }));
	}

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, ndim> coordinates(subkernels.coordinates);
	InitialTemperatureKernel<nodes> temperature(subkernels.temperature);
	IntegrationKernel<nodes, ndim, edim> integration(subkernels.integration);
	ThicknessToNodes<nodes, ndim> thickness(subkernels.thickness);

	struct {
		std::vector<ExternalNodeExpression<ndim, Element>*> node;
		std::vector<ExternalGPsExpression<ndim, Element>*> gp;
	} nonconst;

	for (size_t i = 0; i < subkernels.expressions.node.size(); ++i) {
		ExternalNodeExpression<ndim, Element>* exp = dynamic_cast<ExternalNodeExpression<ndim, Element>*>(subkernels.expressions.node[i]);
		if (subkernels.expressions.node[i]->evaluator->isConst()) {
			for (size_t n = 0; n < nodes; ++n) {
				exp->simd(element, n);
			}
		} else {
			nonconst.node.push_back(exp);
		}
	}

	SIMD volume;
	basis.simd(element);
	for (size_t c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		thickness.simd(element);

		for (size_t i = 0; i < nonconst.node.size(); ++i) {
			for (size_t n = 0; n < nodes; ++n) {
				nonconst.node[i]->simd(element, n);
			}
		}

		temperature.simd(element);

		for (size_t gp = 0; gp < gps; ++gp) {
			integration.simd(element, gp);
			volume = volume + element.det * load1(element.w[gp]);
		}
	}

	subkernels.esize = sizeof(Element);
	for (size_t s = 0; s < SIMD::size; ++s) {
		subkernels.volume += volume[s];
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runElementKernel(const HeatTransferElementOperators &subkernels, SubKernel::Action action)
{
	typedef HeatTransferElement<nodes, gps, ndim, edim> Element; Element element;

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, ndim> coordinates(subkernels.coordinates);
	CoordinatesToGPsKernel<nodes, ndim> coordinatesToGPs(subkernels.coordinates);
	ThicknessToGp<nodes, ndim> thickness(subkernels.thickness);
	TemperatureKernel<nodes> temperature(subkernels.temperature);
	TemperatureToGPsKernel<nodes> temperatureToGPs(subkernels.temperature);
	IntegrationKernel<nodes, ndim, edim> integration(subkernels.integration);
	ConductivityKernel<ndim> conductivity(subkernels.conductivity);
	HeatSourceKernel<nodes> heatSource(subkernels.heatSource);
	AdvectionKernel<nodes, ndim> advection(subkernels.advection);
	MatrixConductivityKernel<nodes, ndim> K(subkernels.K);
	MatrixMassKernel<nodes, 1> M(subkernels.M);
	TemperatureGradientKernel<nodes, gps, ndim> gradient(subkernels.gradient);
	TemperatureFluxKernel<nodes, gps, ndim> flux(subkernels.flux);
	MatricFillerKernel<nodes> outK(subkernels.Kfiller);
	MatricFillerKernel<nodes> outM(subkernels.Mfiller);
	RHSFillerKernel<nodes> outRHS(subkernels.RHSfiller);

	struct {
		std::vector<ExternalNodeExpression<ndim, Element>*> node;
		std::vector<ExternalGPsExpression<ndim, Element>*> gp;
	} nonconst;

	for (size_t i = 0; i < subkernels.expressions.node.size(); ++i) {
		ExternalNodeExpression<ndim, Element>* exp = dynamic_cast<ExternalNodeExpression<ndim, Element>*>(subkernels.expressions.node[i]);
		if (subkernels.expressions.node[i]->evaluator->isConst()) {
			for (size_t n = 0; n < nodes; ++n) {
				exp->simd(element, n);
			}
		} else {
			nonconst.node.push_back(exp);
		}
	}

	for (size_t i = 0; i < subkernels.expressions.gp.size(); ++i) {
		ExternalGPsExpression<ndim, Element>* exp = dynamic_cast<ExternalGPsExpression<ndim, Element>*>(subkernels.expressions.gp[i]);
		if (subkernels.expressions.gp[i]->evaluator->isConst()) {
			for (size_t gp = 0; gp < gps; ++gp) {
				exp->simd(element, gp);
			}
		} else {
			nonconst.gp.push_back(exp);
		}
	}

	// pre-processing of possible constant parameters from ecf
	basis.simd(element);
	conductivity.simd(element);
	thickness.simd(element, 0);

	coordinatesToGPs.setActiveness(action);
	thickness.setActiveness(action);
	temperature.setActiveness(action);
	conductivity.setActiveness(action);
	heatSource.setActiveness(action);
	advection.setActiveness(action);
	K.setActiveness(action);
	M.setActiveness(action);
	gradient.setActiveness(action);
	flux.setActiveness(action);

	outK.setActiveness(action);
	outRHS.setActiveness(action);

	for (size_t c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		if (temperature.isactive) {
			temperature.simd(element);
		}
		for (size_t i = 0; i < nonconst.node.size(); ++i) {
			for (size_t n = 0; n < nodes; ++n) {
				nonconst.node[i]->simd(element, n);
			}
		}

		for (size_t gp = 0; gp < gps; ++gp) {
			integration.simd(element, gp);

			if (coordinatesToGPs.isactive) {
				coordinatesToGPs.simd(element, gp);
			}

			if (temperatureToGPs.isactive) {
				temperatureToGPs.simd(element, gp);
			}

			if (thickness.isactive) {
				thickness.simd(element, gp);
			}

			for (size_t i = 0; i < nonconst.gp.size(); ++i) {
				for (size_t n = 0; n < nodes; ++n) {
					nonconst.gp[i]->simd(element, n);
				}
			}

			if (conductivity.isactive) {
				conductivity.simd(element, gp);
			}
			if (heatSource.isactive) {
				heatSource.simd(element, gp);
			}
			if (advection.isactive) {
				advection.simd(element, gp);
			}
			if (K.isactive) {
				K.simd(element, gp);
			}
			if (M.isactive) {
				M.simd(element, gp);
			}
			if (gradient.isactive) {
				gradient.simd(element, gp);
			}
			if (flux.isactive) {
				flux.simd(element, gp);
			}
		}

		if (gradient.isactive) {
			gradient.store(element);
		}

		if (flux.isactive) {
			flux.store(element);
		}

		if (outK.isactive) {
			outK.simd(element.K);
		}
		if (outM.isactive) {
			outM.simd(element.M);
		}
		if (outRHS.isactive) {
			outRHS.simd(element.f);
		}
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void setBoundaryKernel(HeatTransferBoundaryOperators &subkernels, SubKernel::Action action)
{
	typedef HeatTransferBoundary<nodes, gps, ndim, edim> Element; Element element;

	if (subkernels.heatFlow.expression) {
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.heatFlow.expression->evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.heatFlow[s] = value; }));
	}
	if (subkernels.heatFlux.expression) {
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.heatFlux.expression->evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.heatFlux[s] = value; }));
	}

	if (subkernels.htc.expression) {
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.htc.expression->evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.htc[s] = value; }));
	}
	if (subkernels.externalTemperature.expression) {
		subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
				subkernels.externalTemperature.expression->evaluator,
				[] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.extTemp[s] = value; }));
	}

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, ndim> coordinates(subkernels.coordinates);
	IntegrationKernel<nodes, ndim, edim> integration(subkernels.integration);

	SIMD surface;
	basis.simd(element);
	for (size_t c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		for (size_t gp = 0; gp < gps; ++gp) {
			integration.simd(element, gp);
			surface = surface + element.det * load1(element.w[gp]);
		}
	}

	subkernels.esize = sizeof(Element);
	for (size_t s = 0; s < SIMD::size; ++s) {
		subkernels.surface += surface[s];
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runBoundaryKernel(const HeatTransferBoundaryOperators &subkernels, SubKernel::Action action)
{
	typedef HeatTransferBoundary<nodes, gps, ndim, edim> Element; Element element;

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, ndim> coordinates(subkernels.coordinates);
	CoordinatesToGPsKernel<nodes, ndim> coordinatesToGPs(subkernels.coordinates);
	ThicknessFromNodes<nodes, ndim> thickness(subkernels.thickness);
	ThicknessToGp<nodes, ndim> thicknessToGPs(subkernels.thickness);
	IntegrationKernel<nodes, ndim, edim> integration(subkernels.integration);
	ExternalHeatKernel<nodes> externalHeat(subkernels.externalHeat);
	RHSFillerKernel<nodes> outRHS(subkernels.RHSfiller);

	std::vector<ExternalGPsExpression<ndim, Element>*> nonconst;
	for (size_t i = 0; i < subkernels.expressions.gp.size(); ++i) {
		ExternalGPsExpression<ndim, Element>* exp = dynamic_cast<ExternalGPsExpression<ndim, Element>*>(subkernels.expressions.gp[i]);
		if (subkernels.expressions.gp[i]->evaluator->isConst()) {
			exp->simd(element, 0);
		} else {
			nonconst.push_back(exp);
		}
	}

	basis.simd(element);
	thickness.setActiveness(action);

	outRHS.setActiveness(action);

	for (size_t c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);

		if (thickness.isactive) {
			thickness.simd(element);
		}

		for (size_t gp = 0; gp < gps; ++gp) {
			if (coordinatesToGPs.isactive) {
				coordinatesToGPs.simd(element, gp);
			}

			for (size_t i = 0; i < nonconst.size(); ++i) {
				nonconst[i]->simd(element, gp);
			}

			if (thicknessToGPs.isactive) {
				thicknessToGPs.simd(element, gp);
			}

			integration.simd(element, gp);

			if (externalHeat.isactive) {
				externalHeat.simd(element, gp);
			}
		}

		if (outRHS.isactive) {
			outRHS.simd(element.f);
		}
	}
}

template <size_t ndim>
void setNodeKernel(HeatTransferBoundaryOperators &subkernels, SubKernel::Action action)
{
	typedef HeatTransferNode<ndim> Element; Element element;
	if (subkernels.temperature.expression) {
		auto setter = [] (Element &element, size_t &n, size_t &s, double value) { element.temperature.node[0][s] = element.temperature.initial[0][s] = value; };
		switch (info::mesh->dimension) {
		case 2: subkernels.expressions.node.push_back(new ExternalNodeExpression<2, Element>(subkernels.temperature.expression->evaluator, setter)); break;
		case 3: subkernels.expressions.node.push_back(new ExternalNodeExpression<3, Element>(subkernels.temperature.expression->evaluator, setter)); break;
		}
	}

	CoordinatesKernel<1, ndim> coordinates(subkernels.coordinates);
	InitialTemperatureKernel<1> initTemperature(subkernels.initialTemperature);

	std::vector<ExternalNodeExpression<ndim, Element>*> nonconst;
	for (size_t i = 0; i < subkernels.expressions.node.size(); ++i) {
		ExternalNodeExpression<ndim, Element>* exp = dynamic_cast<ExternalNodeExpression<ndim, Element>*>(subkernels.expressions.node[i]);
		if (subkernels.expressions.node[i]->evaluator->isConst()) {
			exp->simd(element, 0);
		} else {
			nonconst.push_back(exp);
		}
	}

	for (size_t c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		for (size_t i = 0; i < nonconst.size(); ++i) {
			nonconst[i]->simd(element, 0);
		}
		if (initTemperature.isactive) {
			initTemperature.simd(element);
		}
	}
}

template <size_t ndim>
void runNodeKernel(const HeatTransferBoundaryOperators &subkernels, SubKernel::Action action)
{
	typedef HeatTransferNode<ndim> Element; Element element;

	CoordinatesKernel<1, ndim> coordinates(subkernels.coordinates);
	VectorSetterKernel<1, Element> set(subkernels.dirichlet, [] (auto &element, size_t &n, size_t &d, size_t &s) { return element.temperature.node[0][s]; });

	std::vector<ExternalNodeExpression<ndim, Element>*> nonconst;
	for (size_t i = 0; i < subkernels.expressions.node.size(); ++i) {
		ExternalNodeExpression<ndim, Element>* exp = dynamic_cast<ExternalNodeExpression<ndim, Element>*>(subkernels.expressions.node[i]);
		if (subkernels.expressions.node[i]->evaluator->isConst()) {
			exp->simd(element, 0);
		} else {
			nonconst.push_back(exp);
		}
	}

	for (size_t c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		for (size_t i = 0; i < nonconst.size(); ++i) {
			nonconst[i]->simd(element, 0);
		}
		set.simd(element);
	}
}

}



#endif /* SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_KERNEL_H_ */
