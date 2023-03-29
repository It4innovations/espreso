
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

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t ETYPE>
Assembler::measurements StructuralMechanics::manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements,
typename std::enable_if<!
		(ndim == 2 &&
		edim == 2 &&
		ETYPE == StructuralElementType::SYMMETRIC_PLANE), int>::type*
)
{
	return measurements();
}

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t ETYPE>
Assembler::measurements StructuralMechanics::manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements,
typename std::enable_if<
		ndim == 2 &&
		edim == 2 &&
		ETYPE == StructuralElementType::SYMMETRIC_PLANE, int>::type*
)
{
	if (this->K == nullptr) {
		return loop<StructuralMechanicsDataDescriptor, nodes, gps, ndim, edim, ETYPE>(action, ops, elements);
	}
	if ((action == ActionOperator::FILL) || (action == ActionOperator::SOLUTION))
	{
		return conditionsloop<DataDescriptor, nodes, gps, ndim, edim, ETYPE>(action, ops, interval, elements);
	}
	if (elements == 0) return measurements();

	double initStart = eslog::time();

	typename DataDescriptor<nodes, gps, ndim, edim, ETYPE>::Element element;

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if (elements > SIMD::size) {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op)->simd(element);
			} else {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op)->peel(element, elements);
			}
		}
	}

	const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
	bool rotateElasticity = mat->linear_elastic_properties.model != LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC;
	bool constCosSin = mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN;
	// it is dirty hack just to be sure that compiler must assume both variants (currently settings.sigma = 0 and diffusion_split = false)
	bool constElasticity = !settings.contact_interfaces;
	bool constRotation = settings.load_steps == 1;
	if (mat->linear_elastic_properties.model != LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC) {
		if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			if (ndim == 2) {
				rotateElasticity &= mat->coordinate_system.rotation.z.isset;
			}
			if (ndim == 3) {
				rotateElasticity &= mat->coordinate_system.rotation.x.isset | mat->coordinate_system.rotation.y.isset | mat->coordinate_system.rotation.z.isset;
			}
		}
	}

	if (info::ecf->always_update_conductivity) { // TODO
		constElasticity = false;
	}
	bool storeCosSin = settings.reassembling_optimization && action == ActionOperator::ASSEMBLE;
	bool loadCosSin  = settings.reassembling_optimization && action != ActionOperator::ASSEMBLE;

	bool hasAcceleration = false;
	bool constAcceleration = true;
	auto AccelerationEval = configuration.acceleration.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (AccelerationEval != configuration.acceleration.end()) {
		hasAcceleration = true;
		constAcceleration = AccelerationEval->second.x.evaluator != nullptr && AccelerationEval->second.y.evaluator != nullptr && AccelerationEval->second.z.evaluator != nullptr;
	}

	bool hasVelocity = false;
	bool constVelocity = true;
	auto VelocityEval = configuration.angular_velocity.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (VelocityEval != configuration.angular_velocity.end()) {
		hasVelocity = true;
		constVelocity = VelocityEval->second.x.evaluator != nullptr && VelocityEval->second.y.evaluator != nullptr && VelocityEval->second.z.evaluator != nullptr;
	}

	bool constThickness = true;
	auto thicknessEval = settings.thickness.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (thicknessEval != settings.thickness.end()) {
		constThickness = thicknessEval->second.evaluator != nullptr;
	}

	bool cooToGP = mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN || hasVelocity || axisymmetric;
	bool computeK = action == ActionOperator::ASSEMBLE || action == ActionOperator::REASSEMBLE;
	bool computeStress = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.stress;
	bool computeElasticity = computeK | computeStress;
	bool getDisplacement = computeStress;

	esint chunks = elements / SIMD::size;
	double init = eslog::time() - initStart;


	double start = eslog::time();
	switch (action) {
	case ActionOperator::ASSEMBLE  : __SSC_MARK(0xFACE); break;
	case ActionOperator::REASSEMBLE: __SSC_MARK(0xCAFE); break;
	case ActionOperator::SOLUTION  : __SSC_MARK(0xFEED); break; // TODO
	default:
		eslog::error("unsupported action\n");
	}
	// std::cout<<"NDIM  is "<<ndim<< std::endl;
	// std::cout<<"EDIM  is "<<edim<< std::endl;
	// std::cout<<"ETYPE is "<<ETYPE<<std::endl;

	OutputParameterIterator stiffness(this->elements.stiffness, interval);
	OutputParameterIterator rhs(this->elements.rhs, interval);
	auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;

	stiffness += SIMD::size;
	rhs       += SIMD::size;
	for (size_t s = 0; s < SIMD::size; ++s) {
		++procNodes;
	}

	for (esint c = 1; c < chunks; ++c) {

		// coo.simd(element);
		for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t d = 0; d < ndim; ++d) {
					element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
				}
			}
		}
		// integration.simd(element);
		double * __restrict__ out = stiffness.data;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD jacobian0 = zeros(), jacobian1 = zeros(), jacobian2 = zeros(), jacobian3 = zeros();

			for (size_t n = 0; n < nodes; ++n) {
				SIMD coordsX = element.coords[n][0];
				SIMD coordsY = element.coords[n][1];
				SIMD dNX = load1(element.dN[gp][n][0]);
				SIMD dNY = load1(element.dN[gp][n][1]);

				jacobian0 = jacobian0 + dNX * coordsX;
				jacobian1 = jacobian1 + dNX * coordsY;
				jacobian2 = jacobian2 + dNY * coordsX;
				jacobian3 = jacobian3 + dNY * coordsY;
			}

			element.det[gp] = jacobian0 * jacobian3 - jacobian1 * jacobian2;

			SIMD detJx = ones() / element.det[gp];
			SIMD inv0 =  detJx * jacobian3;
			SIMD inv1 = -detJx * jacobian1;
			SIMD inv2 = -detJx * jacobian2;
			SIMD inv3 =  detJx * jacobian0;

			for (size_t n = 0; n < nodes; ++n) {
				SIMD dNX = load1(element.dN[gp][n][0]);
				SIMD dNY = load1(element.dN[gp][n][1]);
				element.dND[gp][n][0] = inv0 * dNX + inv1 * dNY;
				element.dND[gp][n][1] = inv2 * dNX + inv3 * dNY;
			}

			// stiffness.simd(element);
			SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]);
			SIMD c00 = element.elasticity[gp][0];
			SIMD c01 = element.elasticity[gp][1], c11 = element.elasticity[gp][3];
			SIMD c02 = element.elasticity[gp][2], c12 = element.elasticity[gp][4], c22 = element.elasticity[gp][5];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD a = nx * c00 + ny * c02;
					SIMD c = nx * c02 + ny * c22;

					size_t i = n * (2 * nodes) + m - ((n + 1) * n / 2);
					SIMD xx = load(out + i * SIMD::size);
					xx = xx + scale * (a * mx + c * my);
					store(out + i * SIMD::size, xx);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD b = nx * c01 + ny * c12;
					SIMD c = nx * c02 + ny * c22;

					size_t i =  n * (2 * nodes) + (m + nodes) - ((n + 1) * n / 2);
					SIMD xy = load(out + i * SIMD::size);
					xy = xy + scale * (b * my + c * mx);
					store(out + i * SIMD::size, xy);
				}
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD b = ny * c11 + nx * c12;
					SIMD c = ny * c12 + nx * c22;

					size_t i = (n + nodes) * (2 * nodes) + (m + nodes) - (((n + nodes) + 1) * (n + nodes) / 2);
					SIMD yy = load(out + i * SIMD::size);
					yy = yy + scale * (b * my + c * mx);
					store(out + i * SIMD::size, yy);
				}
			}
		}
		// acceleration.simd(element);
		double * __restrict__ outrhs = rhs.data;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(outrhs + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(outrhs + (1 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]) * element.ecf.density[gp] * load1(element.N[gp][n]);
				fx = fx + scale * element.ecf.acceleration[gp][0];
				fy = fy + scale * element.ecf.acceleration[gp][1];
			}
			store(outrhs + (0 * nodes + n) * SIMD::size, fx);
			store(outrhs + (1 * nodes + n) * SIMD::size, fy);
		}

		stiffness += SIMD::size;
		rhs       += SIMD::size;
	}

	switch (action) {
	case ActionOperator::ASSEMBLE  : __SSC_MARK(0xDEAD); break;
	case ActionOperator::REASSEMBLE: __SSC_MARK(0xFADE); break;
	case ActionOperator::SOLUTION  : __SSC_MARK(0xBEED); break; // TODO
	default:
		eslog::error("unsupported action\n");
	}
	double loop = eslog::time() - start;

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
	return measurements(init, loop);
}

template <int etype>
Assembler::measurements StructuralMechanics::instantiateManual2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return manualloop<StructuralMechanicsDataDescriptor, 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return manualloop<StructuralMechanicsDataDescriptor, 6, StructuralMechanicsGPC::TRIANGLE6, 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return manualloop<StructuralMechanicsDataDescriptor, 4, StructuralMechanicsGPC::SQUARE4  , 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return manualloop<StructuralMechanicsDataDescriptor, 8, StructuralMechanicsGPC::SQUARE8  , 2, 2, etype>(action, ops, interval, elements); break;
	default: return measurements();
	}
}

template <int etype>
Assembler::measurements StructuralMechanics::instantiateManual3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TETRA4):    return manualloop<StructuralMechanicsDataDescriptor,  4, StructuralMechanicsGPC::TETRA4    , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::TETRA10):   return manualloop<StructuralMechanicsDataDescriptor, 10, StructuralMechanicsGPC::TETRA10   , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return manualloop<StructuralMechanicsDataDescriptor,  5, StructuralMechanicsGPC::PYRAMID5  , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): return manualloop<StructuralMechanicsDataDescriptor, 13, StructuralMechanicsGPC::PYRAMID13 , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA6):   return manualloop<StructuralMechanicsDataDescriptor,  6, StructuralMechanicsGPC::PRISMA6   , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA15):  return manualloop<StructuralMechanicsDataDescriptor, 15, StructuralMechanicsGPC::PRISMA15  , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA8):     return manualloop<StructuralMechanicsDataDescriptor,  8, StructuralMechanicsGPC::HEXA8     , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA20):    return manualloop<StructuralMechanicsDataDescriptor, 20, StructuralMechanicsGPC::HEXA20    , 3, 3, etype>(action, ops, interval, elements); break;
	default: return measurements();
	}
}

Assembler::measurements StructuralMechanics::instantiateManual(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (etype) {
		// elements
		case StructuralMechanicsElementType::SYMMETRIC_PLANE             : return instantiateManual2D<StructuralMechanicsElementType::SYMMETRIC_PLANE             >(action, code, ops, interval, elements);
		case StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC: return instantiateManual2D<StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC>(action, code, ops, interval, elements);
		}
	case 3:
		switch (etype) {
		// elements
		case StructuralMechanicsElementType::SYMMETRIC_VOLUME: return instantiateManual3D<StructuralMechanicsElementType::SYMMETRIC_VOLUME>(action, code, ops, interval, elements);
		}
	}
	return measurements();
}

}

