
#include "structuralmechanics.h"
#include "structuralmechanics.generator.h"
#include "assembler.hpp"

#include "analysis/assembler/operators/info.h"
#include "analysis/assembler/operators/basis.h"
#include "analysis/assembler/operators/coordinates.h"
#include "analysis/assembler/operators/expression.h"
#include "analysis/assembler/operators/integration.h"
#include "analysis/assembler/operators/structuralmechanics.f.h"
#include "analysis/assembler/operators/structuralmechanics.K.h"
#include "analysis/assembler/operators/filler.h"

#include "basis/expression/variable.h"
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

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct updateAcceleration {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element, Evaluator* evaluator)
	{
		double results[SIMD::size * gps];
		evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.ecf.acceleration[gp][0][s] = results[gps * s + gp];
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype> struct updateVelocity;

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct updateVelocity<DataDescriptor, nodes, gps, 2, edim, etype> {
	void operator()(typename DataDescriptor<nodes, gps, 2, edim, etype>::Element &element, Evaluator* evaluator)
	{
		double results[SIMD::size * gps];
		evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.ecf.angularVelocity[gp][s] = results[gps * s + gp];
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct updateVelocity<DataDescriptor, nodes, gps, 3, edim, etype> {
	void operator()(typename DataDescriptor<nodes, gps, 3, edim, etype>::Element &element, Evaluator* evaluator)
	{
		double results[SIMD::size * gps];
		evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.ecf.angularVelocity[gp][0][s] = results[gps * s + gp];
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
Assembler::measurements StructuralMechanics::conditionsloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	double initStart, initEnd;
	initStart = eslog::time();
	if (this->K == nullptr) {
		return loop<StructuralMechanicsDataDescriptor, nodes, gps, ndim, edim, etype>(action, ops, elements);
	}
	if (elements == 0) return { .0, .0 };

	typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element element;

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if (elements > SIMD::size) {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, etype>*>(*op)->simd(element);
			} else {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, etype>*>(*op)->peel(element, elements);
			}
		}
	}

	auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;
	CoordinatesToElementNodes<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > coo(interval, procNodes);
	CoordinatesToElementNodesAndGPs<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > cooAndGps(interval, procNodes);
	Integration<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > integration(interval);
	StructuralMechanicsStiffness<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > stiffness(interval, this->elements.stiffness);
	Acceleration<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > acceleration(interval, this->elements.rhs);
	AngularVelocity<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > velocity(interval, this->elements.rhs);

//	applyRotation<DataDescriptor, nodes, gps, ndim, edim, etype> rotation(interval);
	updateAcceleration<DataDescriptor, nodes, gps, ndim, edim, etype> updateAcc;
	updateVelocity<DataDescriptor, nodes, gps, ndim, edim, etype> updateVelocity;

	SymmetricMatricFiller<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > upperFiller(interval, ndim, this->elements.stiffness, this->K);
	GeneralMatricFiller<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > fullFiller(interval, ndim, this->elements.stiffness, this->K);
	VectorFiller<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > rhsFiller(interval, ndim, this->elements.rhs, this->f);

//	TemperatureGradient<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > gradient(interval, Results::gradient);
//	TemperatureFlux<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > flux(interval, Results::flux);

	coo.move(SIMD::size);
	cooAndGps.move(SIMD::size);
//	temp.move(SIMD::size);
	stiffness.move(SIMD::size);
	acceleration.move(SIMD::size);
	velocity.move(SIMD::size);
	upperFiller.move(SIMD::size);
	fullFiller.move(SIMD::size);
	rhsFiller.move(SIMD::size);
//	heatSource.move(SIMD::size);

//	gradient.move(ndim * SIMD::size);
//	flux.move(ndim * SIMD::size);

	const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
	bool rotateConductivity = mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	// it is dirty hack just to be sure that compiler must assume both variants (currently settings.sigma = 0 and diffusion_split = false)
//	bool constConductivity = !settings.diffusion_split;
//	bool constRotation = settings.sigma == 0;
//	if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
//		if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
//			if (ndim == 2) {
//				rotateConductivity &= mat->coordinate_system.rotation.z.isset;
//			}
//			if (ndim == 3) {
//				rotateConductivity &= mat->coordinate_system.rotation.x.isset | mat->coordinate_system.rotation.y.isset | mat->coordinate_system.rotation.z.isset;
//			}
//		}
//	}

	bool hasAcceleration = false;
	bool constAcceleration = true;
	auto AccelerationEval = configuration.acceleration.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (AccelerationEval != configuration.acceleration.end()) {
		hasAcceleration = true;
		constAcceleration = AccelerationEval->second.x.evaluator->params.general.size() == 0 && AccelerationEval->second.y.evaluator->params.general.size() == 0 && AccelerationEval->second.z.evaluator->params.general.size() == 0;
	}

	bool hasVelocity = false;
	bool constVelocity = true;
	auto VelocityEval = configuration.acceleration.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (AccelerationEval != configuration.acceleration.end()) {
		hasAcceleration = true;
		constAcceleration = AccelerationEval->second.x.evaluator->params.general.size() == 0 && AccelerationEval->second.y.evaluator->params.general.size() == 0 && AccelerationEval->second.z.evaluator->params.general.size() == 0;
	}

	bool cooToGP = mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN;
	bool computeK = action == ActionOperator::ASSEMBLE || action == ActionOperator::REASSEMBLE;
	bool computeGradient = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.gradient;
	bool computeFlux = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.flux;
	bool computeConductivity = computeK | computeFlux;
	bool getTemp = computeGradient || computeFlux;
	bool isfullMatrix = this->K->shape == Matrix_Shape::FULL;

	initEnd = eslog::time();
	double start, end;

	if (action == ActionOperator::Action::ASSEMBLE)
	{
		start = eslog::time();
		__SSC_MARK(0xFACE);
		esint chunks = elements / SIMD::size;
		for (esint c = 1; c < chunks; ++c) {

		}
		__SSC_MARK(0xDEAD);
		end = eslog::time();
	}

	if (action == ActionOperator::Action::REASSEMBLE)
	{
		start = eslog::time();
		__SSC_MARK(0xCAFE);
		esint chunks = elements / SIMD::size;
		for (esint c = 1; c < chunks; ++c) {

		}
		__SSC_MARK(0xDADE);
		end = eslog::time();
	}

	if (action != ActionOperator::Action::REASSEMBLE && action != ActionOperator::Action::ASSEMBLE)
	{
		start = eslog::time();
		__SSC_MARK(0xCAFE);
		esint chunks = elements / SIMD::size;
		for (esint c = 1; c < chunks; ++c) {

		}
		__SSC_MARK(0xDADE);
		end = eslog::time();
	}

	esint chunks = elements / SIMD::size;
	for (esint c = 1; c < chunks; ++c) {
		if (cooToGP) {
			cooAndGps.simd(element);
		} else {
			coo.simd(element);
		}
		integration.simd(element);
//		if (getTemp) {
//			temp.simd(element);
//		}
		if (computeConductivity) {
//			if (!constConductivity) {
//				updateConductivity<DataDescriptor, nodes, gps, ndim, edim, etype>()(element, mat);
//			}
//			if (rotateConductivity) {
//				if (!constRotation) {
//					updateRotation<DataDescriptor, nodes, gps, ndim, edim, etype>()(element, mat);
//				}
//				rotation(element, mat);
//			}
		}

		if (computeK) {
			stiffness.simd(element);
			if (hasAcceleration) {
				if (!constAcceleration) {
					updateAcc(element, AccelerationEval->second.x.evaluator);
				}
				acceleration.simd(element);
			}
		}
		if (action == ActionOperator::FILL) {
			if (isfullMatrix) {
				fullFiller.simd(element);
			} else {
				upperFiller.simd(element);
			}
			rhsFiller.simd(element);
		}
	}

	if (elements % SIMD::size) {
		eslog::error("peel loop is not supported\n");
		// peel is never needed
	}

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if ((*op)->isconst) {
				(*op)->move(-(int)std::min(elements, (esint)SIMD::size));
			} else {
				(*op)->move(-(esint)SIMD::size);
			}
		}
	}
	return {initEnd - initStart, end - start};
}

template <int etype>
Assembler::measurements StructuralMechanics::instantiateConditions2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return conditionsloop<StructuralMechanicsDataDescriptor, 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return conditionsloop<StructuralMechanicsDataDescriptor, 6, StructuralMechanicsGPC::TRIANGLE6, 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return conditionsloop<StructuralMechanicsDataDescriptor, 4, StructuralMechanicsGPC::SQUARE4  , 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return conditionsloop<StructuralMechanicsDataDescriptor, 8, StructuralMechanicsGPC::SQUARE8  , 2, 2, etype>(action, ops, interval, elements); break;
	default: return {0.0, 0.0};
	}
}

template <int etype>
Assembler::measurements StructuralMechanics::instantiateConditions3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TETRA4):    return conditionsloop<StructuralMechanicsDataDescriptor,  4, StructuralMechanicsGPC::TETRA4    , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::TETRA10):   return conditionsloop<StructuralMechanicsDataDescriptor, 10, StructuralMechanicsGPC::TETRA10   , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return conditionsloop<StructuralMechanicsDataDescriptor,  5, StructuralMechanicsGPC::PYRAMID5  , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): return conditionsloop<StructuralMechanicsDataDescriptor, 13, StructuralMechanicsGPC::PYRAMID13 , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA6):   return conditionsloop<StructuralMechanicsDataDescriptor,  6, StructuralMechanicsGPC::PRISMA6   , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA15):  return conditionsloop<StructuralMechanicsDataDescriptor, 15, StructuralMechanicsGPC::PRISMA15  , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA8):     return conditionsloop<StructuralMechanicsDataDescriptor,  8, StructuralMechanicsGPC::HEXA8     , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA20):    return conditionsloop<StructuralMechanicsDataDescriptor, 20, StructuralMechanicsGPC::HEXA20    , 3, 3, etype>(action, ops, interval, elements); break;
	default: return {0.0, 0.0};
	}
}

Assembler::measurements StructuralMechanics::instantiateConditions(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (etype) {
		// elements
		case StructuralMechanicsElementType::SYMMETRIC_PLANE             : return instantiateConditions2D<StructuralMechanicsElementType::SYMMETRIC_PLANE             >(action, code, ops, interval, elements);
		case StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC: return instantiateConditions2D<StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC>(action, code, ops, interval, elements);
		}
	case 3:
		switch (etype) {
		// elements
		case StructuralMechanicsElementType::SYMMETRIC_VOLUME: return instantiateConditions3D<StructuralMechanicsElementType::SYMMETRIC_VOLUME>(action, code, ops, interval, elements);
		}
	}
	return {0.0, 0.0};
}

}

