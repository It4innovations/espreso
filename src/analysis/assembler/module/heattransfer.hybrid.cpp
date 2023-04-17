
#include "heattransfer.h"
#include "heattransfer.element.h"
#include "analysis/assembler/module/assembler.hpp"

#include "analysis/assembler/subkernel/basis.h"
#include "analysis/assembler/subkernel/coordinates.h"
#include "analysis/assembler/subkernel/temperature.h"
#include "analysis/assembler/subkernel/integration.h"
#include "analysis/assembler/subkernel/expression.h"
#include "analysis/assembler/subkernel/heattransfer/advection.h"
#include "analysis/assembler/subkernel/heattransfer/conductivity.h"
#include "analysis/assembler/subkernel/heattransfer/coordinatesystem.h"
#include "analysis/assembler/subkernel/heattransfer/flux.h"
#include "analysis/assembler/subkernel/heattransfer/gradient.h"
#include "analysis/assembler/subkernel/heattransfer/matrix.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

#include "analysis/scheme/steadystate.h"
#include "math/physics/matrix_distributed.h"

#include <numeric>
#include <algorithm>

namespace espreso {

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype> struct HeatTransferUserParameters {
	static void analyze(HeatTransfer::SubKernels &subkernels) {}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct HeatTransferUserParameters<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC> {

	static void analyze(HeatTransfer::SubKernels &subkernels)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct HeatTransferUserParameters<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {

	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		typedef DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> Physics;

		if (subkernels.translation.cooSystem) {
			switch (subkernels.translation.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->rotation.z.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->center.x.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->center.y.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][1][s] = value; }));
				break;
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct HeatTransferUserParameters<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC> {

	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		typedef DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC> Physics;

		if (subkernels.advection.expression) {
			subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
					subkernels.advection.expression->x.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
					subkernels.advection.expression->y.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][1][s] = value; }));
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct HeatTransferUserParameters<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {

	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		typedef DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> Physics;

		if (subkernels.translation.cooSystem) {
			if (subkernels.translation.cooSystem) {
				switch (subkernels.translation.type) {
				case CoordinateSystemConfiguration::TYPE::CARTESIAN:
					subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
							subkernels.translation.cooSystem->rotation.z.evaluator,
							[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
					break;
				case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
					subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
							subkernels.translation.cooSystem->center.x.evaluator,
							[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
					subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
							subkernels.translation.cooSystem->center.y.evaluator,
							[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][1][s] = value; }));
					break;
				}
			}
		}
		if (subkernels.advection.expression) {
			subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
					subkernels.advection.expression->x.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
					subkernels.advection.expression->y.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][1][s] = value; }));
		}
	}
};


template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct HeatTransferUserParameters<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC> {

	static void analyze(HeatTransfer::SubKernels &subkernels)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct HeatTransferUserParameters<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {

	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		typedef DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> Physics;

		if (subkernels.translation.cooSystem) {
			switch (subkernels.translation.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->rotation.x.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->rotation.y.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][1][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->rotation.z.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][2][s] = value; }));
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			case CoordinateSystemConfiguration::TYPE::SPHERICAL:
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->center.x.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->center.y.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][1][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->center.z.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][2][s] = value; }));
				break;
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct HeatTransferUserParameters<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC> {

	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		typedef DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC> Physics;

		if (subkernels.advection.expression) {
			subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
					subkernels.advection.expression->x.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
					subkernels.advection.expression->y.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
					subkernels.advection.expression->z.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][2][s] = value; }));
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct HeatTransferUserParameters<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {

	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		typedef DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> Physics;

		if (subkernels.translation.cooSystem) {
			switch (subkernels.translation.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->rotation.x.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->rotation.y.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][1][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->rotation.z.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][2][s] = value; }));
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			case CoordinateSystemConfiguration::TYPE::SPHERICAL:
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->center.x.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->center.y.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][1][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
						subkernels.translation.cooSystem->center.z.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][2][s] = value; }));
				break;
			}
		}
		if (subkernels.advection.expression) {
			subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
					subkernels.advection.expression->x.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
					subkernels.advection.expression->y.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpressionKernel<gps, Physics>(
					subkernels.advection.expression->z.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][2][s] = value; }));
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
Assembler::measurements HeatTransfer::hybridloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	if (elements == 0) return measurements();

	if (action == ActionOperator::PREPROCESS) {
		HeatTransferUserParameters<DataDescriptor, nodes, gps, ndim, edim, etype>::analyze(subkernels[interval]);
		return measurements();
	}

	double initStart = eslog::time();

	typedef DataDescriptor<nodes, gps, ndim, edim, etype> Physics;
	typename Physics::Element element;

	BasisElementKernel<code, nodes, gps, edim, Physics> basis(subkernels[interval].basis);
	CoordinatesKernelToElementNodes<nodes, gps, ndim, Physics> coordinates(subkernels[interval].coordinates);
	TemperatureSubKernelToElementNodes<nodes, gps, Physics> temperature(subkernels[interval].temperature);
	IntegrationSubKernel<nodes, gps, ndim, edim, Physics> integration(subkernels[interval].integration);
	ConductivitySubKernel<gps, ndim, etype, Physics> conductivity(subkernels[interval].conductivity);
	HeatTransferCoordinateSystemSubKernel<gps, ndim, etype, Physics> translation(subkernels[interval].translation);
	AdvectionSubkernel<nodes, gps, ndim, etype, Physics> advection(subkernels[interval].advection);
	HeatTransferStiffnessKernel<nodes, gps, ndim, etype, Physics> K(subkernels[interval].K);
	TemperatureGradient<nodes, gps, ndim, Physics> gradient(subkernels[interval].gradient);
	TemperatureFlux<nodes, gps, ndim, etype, Physics> flux(subkernels[interval].flux);

	std::vector<ExternalGPsExpressionKernel<gps, Physics>*> nonconst;
	for (size_t i = 0; i < subkernels[interval].expressions.size(); ++i) {
		if (subkernels[interval].expressions[i]->evaluator->isConst()) {
			dynamic_cast<ExternalGPsExpressionKernel<gps, Physics>*>(subkernels[interval].expressions[i])->simd(element);
		} else {
			nonconst.push_back(dynamic_cast<ExternalGPsExpressionKernel<gps, Physics>*>(subkernels[interval].expressions[i]));
		}
	}

	basis.simd(element);
	conductivity.simd(element);
	if (translation.isactive) {
		translation.simd(element);
	}

	conductivity.setActiveness(!conductivity.isconst, action);
	translation.setActiveness(!translation.isconst, action);
	advection.setActiveness(action);
	temperature.setActiveness(action);
	K.setActiveness(action);
	gradient.setActiveness(action);
	flux.setActiveness(action);

	esint chunks = (elements - 1) / SIMD::size + 1;
	double init = eslog::time() - initStart;

	double start = eslog::time();
	switch (action) {
	case ActionOperator::ASSEMBLE  : __SSC_MARK(0xFACE); break;
	case ActionOperator::REASSEMBLE: __SSC_MARK(0xCAFE); break;
	case ActionOperator::SOLUTION  : __SSC_MARK(0xFEED); break;
	default:
		eslog::error("unsupported action\n");
	}

//	printf("sub-kernels: ");
	for (esint c = 0; c < chunks; ++c) {
		coordinates.simd(element);
//		if (c == 0) printf("coordinates ");
		if (temperature.isactive) {
			temperature.simd(element);
//			if (c == 0) printf("temp ");
		}
		integration.simd(element);
//		if (c == 0) printf("integrate ");
		if (conductivity.isactive) {
			conductivity.simd(element);
//			if (c == 0) printf("conductivity ");
		}
		if (translation.isactive) {
			translation.simd(element);
//			if (c == 0) printf("translation ");
		}
		if (advection.isactive) {
			advection.simd(element);
		}
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

	switch (action) {
	case ActionOperator::ASSEMBLE  : __SSC_MARK(0xDEAD); break;
	case ActionOperator::REASSEMBLE: __SSC_MARK(0xFADE); break;
	case ActionOperator::SOLUTION  : __SSC_MARK(0xBEED); break;
	default:
		eslog::error("unsupported action\n");
	}
	double loop = eslog::time() - start;

	return measurements(init, loop);
}

template <int etype>
Assembler::measurements HeatTransfer::instantiateHybrid2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return hybridloop<HeatTransferDataDescriptor, Element::CODE::TRIANGLE3, 3, HeatTransferGPC::TRIANGLE3, 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return hybridloop<HeatTransferDataDescriptor, Element::CODE::TRIANGLE6, 6, HeatTransferGPC::TRIANGLE6, 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE4  ): return hybridloop<HeatTransferDataDescriptor, Element::CODE::SQUARE4  , 4, HeatTransferGPC::SQUARE4  , 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE8  ): return hybridloop<HeatTransferDataDescriptor, Element::CODE::SQUARE8  , 8, HeatTransferGPC::SQUARE8  , 2, 2, etype>(action, ops, interval, elements); break;
	default: return measurements();
	};
}

template <int etype>
Assembler::measurements HeatTransfer::instantiateHybrid3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TETRA4   ): return hybridloop<HeatTransferDataDescriptor, Element::CODE::TETRA4   ,  4, HeatTransferGPC::TETRA4    , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::TETRA10  ): return hybridloop<HeatTransferDataDescriptor, Element::CODE::TETRA10  , 10, HeatTransferGPC::TETRA10   , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5 ): return hybridloop<HeatTransferDataDescriptor, Element::CODE::PYRAMID5 ,  5, HeatTransferGPC::PYRAMID5  , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): return hybridloop<HeatTransferDataDescriptor, Element::CODE::PYRAMID13, 13, HeatTransferGPC::PYRAMID13 , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA6  ): return hybridloop<HeatTransferDataDescriptor, Element::CODE::PRISMA6  ,  6, HeatTransferGPC::PRISMA6   , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA15 ): return hybridloop<HeatTransferDataDescriptor, Element::CODE::PRISMA15 , 15, HeatTransferGPC::PRISMA15  , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA8    ): return hybridloop<HeatTransferDataDescriptor, Element::CODE::HEXA8    ,  8, HeatTransferGPC::HEXA8     , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA20   ): return hybridloop<HeatTransferDataDescriptor, Element::CODE::HEXA20   , 20, HeatTransferGPC::HEXA20    , 3, 3, etype>(action, ops, interval, elements); break;
	default: return measurements();
	};
}

Assembler::measurements HeatTransfer::instantiateHybrid(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (etype) {
		// elements
		case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return instantiateHybrid2D<HeatTransferElementType::SYMMETRIC_ISOTROPIC >(action, code, ops, interval, elements);
		case HeatTransferElementType::SYMMETRIC_GENERAL   : return instantiateHybrid2D<HeatTransferElementType::SYMMETRIC_GENERAL   >(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return instantiateHybrid2D<HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_GENERAL  : return instantiateHybrid2D<HeatTransferElementType::ASYMMETRIC_GENERAL  >(action, code, ops, interval, elements);
		}
	case 3:
		switch (etype) {
		// elements
		case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return instantiateHybrid3D<HeatTransferElementType::SYMMETRIC_ISOTROPIC >(action, code, ops, interval, elements);
		case HeatTransferElementType::SYMMETRIC_GENERAL   : return instantiateHybrid3D<HeatTransferElementType::SYMMETRIC_GENERAL   >(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return instantiateHybrid3D<HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_GENERAL  : return instantiateHybrid3D<HeatTransferElementType::ASYMMETRIC_GENERAL  >(action, code, ops, interval, elements);
		}
	}
	return measurements();
}

}
