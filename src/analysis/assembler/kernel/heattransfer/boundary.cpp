
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/assembler/module/assembler.hpp"
#include "analysis/assembler/subkernel/heattransfer/externalheat.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

#include "math/physics/matrix_distributed.h"

#include <numeric>
#include <algorithm>

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void preprocess(HeatTransferBoundarySubKernelsList &subkernels)
{
	typedef HeatTransferBoundaryDescriptor<nodes, gps, ndim, edim> Physics;
	typename Physics::Element element;
	if (subkernels.heatFlow.expression) {
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.heatFlow.expression->evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.heatFlow[gp][s] = value; }));
	}
	if (subkernels.heatFlux.expression) {
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.heatFlux.expression->evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.heatFlux[gp][s] = value; }));
	}

	if (subkernels.htc.expression) {
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.htc.expression->evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.htc[gp][s] = value; }));
	}
	if (subkernels.externalTemperature.expression) {
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.externalTemperature.expression->evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.extTemp[gp][s] = value; }));
	}

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, gps, ndim, Physics> coordinates(subkernels.coordinates);
	IntegrationKernel<nodes, gps, ndim, edim, Physics> integration(subkernels.integration);

	basis.simd(element);
	SIMD surface;
	for (esint c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		integration.simd(element);
		for (size_t gp = 0; gp < gps; ++gp) {
			surface = surface + element.det[gp] * load1(element.w[gp]);
		}
	}

	subkernels.esize = sizeof(typename Physics::Element);
	for (size_t s = 0; s < SIMD::size; ++s) {
		subkernels.surface += surface[s];
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void compute(const HeatTransferBoundarySubKernelsList &subkernels, Assembler::Action action)
{
	typedef HeatTransferBoundaryDescriptor<nodes, gps, ndim, edim> Physics;
	typename Physics::Element element;

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, gps, ndim, Physics> coordinates(subkernels.coordinates);
	ThicknessFromNodes<nodes, ndim, Physics> thickness(subkernels.thickness);
	IntegrationKernel<nodes, gps, ndim, edim, Physics> integration(subkernels.integration);
	ExternalHeatKernel<nodes, gps, ndim, Physics> externalHeat(subkernels.externalHeat);

	std::vector<ExternalGPsExpression<gps, Physics>*> nonconst;
	for (size_t i = 0; i < subkernels.expressions.size(); ++i) {
		if (subkernels.expressions[i]->evaluator->isConst()) {
			dynamic_cast<ExternalGPsExpression<gps, Physics>*>(subkernels.expressions[i])->simd(element);
		} else {
			nonconst.push_back(dynamic_cast<ExternalGPsExpression<gps, Physics>*>(subkernels.expressions[i]));
		}
	}

	basis.simd(element);
	thickness.setActiveness(action);

	for (esint c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
//		if (c == 0) printf("coordinates ");
		if (thickness.isactive) {
			thickness.simd(element);
//			if (c == 0) printf("thickness ");
		}
		integration.simd(element);
//		if (c == 0) printf("integration ");
		if (externalHeat.isactive) {
			externalHeat.simd(element);
//			if (c == 0) printf("heat ");
		}
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void fill(const HeatTransferBoundarySubKernelsList &subkernels)
{
	typedef HeatTransferBoundaryDescriptor<nodes, gps, ndim, edim> Physics;
	typename Physics::Element element;

	VectorFillerKernel<nodes, Physics> RHS(subkernels.RHSfiller);

	for (esint c = 0; c < subkernels.chunks; ++c) {
		if (RHS.isactive) {
			RHS.simd(element);
		}
	}
}

template <size_t ndim>
void initDirichlet(HeatTransferBoundarySubKernelsList &subkernels)
{
	typedef HeatTransferBoundaryDescriptor<1, 1, ndim, 0> Physics;
	if (subkernels.temperature.expression) {
		subkernels.expressions.push_back(new ExternalNodeExpression<1, Physics>(
				subkernels.temperature.expression->evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.temp[0][s] = value; }));
	}
}

template <size_t ndim>
void dirichlet(const HeatTransferBoundarySubKernelsList &subkernels)
{
	typedef HeatTransferBoundaryDescriptor<1, 1, ndim, 0> Physics;
	typename Physics::Element element;

	CoordinatesKernel<1, 1, ndim, Physics> coordinates(subkernels.coordinates);
	VectorSetterKernel<1, Physics> set(subkernels.dirichlet, [] (auto &element, size_t &n, size_t &d, size_t &s) { return element.temp[n][s]; });

	std::vector<ExternalNodeExpression<1, Physics>*> nonconst;
	for (size_t i = 0; i < subkernels.expressions.size(); ++i) {
		if (subkernels.expressions[i]->evaluator->isConst()) {
			dynamic_cast<ExternalNodeExpression<1, Physics>*>(subkernels.expressions[i])->simd(element);
		} else {
			nonconst.push_back(dynamic_cast<ExternalNodeExpression<1, Physics>*>(subkernels.expressions[i]));
		}
	}

	for (esint c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		for (size_t i = 0; i < nonconst.size(); ++i) {
			nonconst[i]->simd(element);
		}
		set.simd(element);
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runAction(HeatTransferBoundarySubKernelsList &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: preprocess<code, nodes, gps, ndim, edim>(subkernels); break;
	case Assembler::Action::FILL: fill<code, nodes, gps, ndim, edim>(subkernels); break;
	default: compute<code, nodes, gps, ndim, edim>(subkernels, action); break;
	}
}

template <size_t ndim, size_t edim> void addBC(HeatTransferBoundarySubKernelsList &subkernels, Assembler::Action action);

template <> void addBC<2, 1>(HeatTransferBoundarySubKernelsList &subkernels, Assembler::Action action)
{
	switch (subkernels.code) {
	case static_cast<size_t>(Element::CODE::LINE2): runAction<Element::CODE::LINE2, 2, HeatTransferGPC::LINE2, 2, 1>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::LINE3): runAction<Element::CODE::LINE3, 3, HeatTransferGPC::LINE3, 2, 1>(subkernels, action); break;
	}
}

template <> void addBC<2, 0>(HeatTransferBoundarySubKernelsList &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: initDirichlet<2>(subkernels); break;
	case Assembler::Action::ASSEMBLE: case Assembler::Action::REASSEMBLE: dirichlet<2>(subkernels); break;
	}
}
template <> void addBC<3, 2>(HeatTransferBoundarySubKernelsList &subkernels, Assembler::Action action)
{
	switch (subkernels.code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): runAction<Element::CODE::TRIANGLE3, 3, HeatTransferGPC::TRIANGLE3, 3, 2>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): runAction<Element::CODE::TRIANGLE6, 6, HeatTransferGPC::TRIANGLE6, 3, 2>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::SQUARE4  ): runAction<Element::CODE::SQUARE4  , 4, HeatTransferGPC::SQUARE4  , 3, 2>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::SQUARE8  ): runAction<Element::CODE::SQUARE8  , 8, HeatTransferGPC::SQUARE8  , 3, 2>(subkernels, action); break;
	}
}

template <> void addBC<3, 1>(HeatTransferBoundarySubKernelsList &subkernels, Assembler::Action action)
{
	switch (subkernels.code) {
	case static_cast<size_t>(Element::CODE::LINE2): runAction<Element::CODE::LINE2, 2, HeatTransferGPC::LINE2, 3, 1>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::LINE3): runAction<Element::CODE::LINE3, 3, HeatTransferGPC::LINE3, 3, 1>(subkernels, action); break;
	}
}

template <> void addBC<3, 0>(HeatTransferBoundarySubKernelsList &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: initDirichlet<3>(subkernels); break;
	case Assembler::Action::ASSEMBLE: case Assembler::Action::REASSEMBLE: dirichlet<3>(subkernels); break;
	}
}

void HeatTransfer::runBoundary(Action action, size_t region, size_t interval)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (info::mesh->boundaryRegions[region]->dimension) {
		case 0: addBC<2, 0>(boundary[region][interval], action); break;
		case 1: addBC<2, 1>(boundary[region][interval], action); break;
		} break;
	case 3:
		switch (info::mesh->boundaryRegions[region]->dimension) {
		case 0: addBC<3, 0>(boundary[region][interval], action); break;
		case 1: addBC<3, 1>(boundary[region][interval], action); break;
		case 2: addBC<3, 2>(boundary[region][interval], action); break;
		} break;
	}
}

}


