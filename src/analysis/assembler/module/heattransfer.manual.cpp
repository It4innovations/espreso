
#include "heattransfer.h"
#include "heattransfer.generator.h"
#include "assembler.hpp"

#include "analysis/assembler/operators/info.h"
#include "analysis/assembler/operators/basis.h"
#include "analysis/assembler/operators/coordinates.h"
#include "analysis/assembler/operators/temperature.h"
#include "analysis/assembler/operators/advection.h"
#include "analysis/assembler/operators/expression.h"
#include "analysis/assembler/operators/integration.h"
#include "analysis/assembler/operators/conductivity.coordinatesystem.h"
#include "analysis/assembler/operators/heattransfer.f.h"
#include "analysis/assembler/operators/heattransfer.K.h"
#include "analysis/assembler/operators/filler.h"
#include "analysis/assembler/operators/gradient.h"
#include "analysis/assembler/operators/flux.h"

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
#include <iostream>

namespace espreso {


template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t ETYPE>
Assembler::measurements HeatTransfer::manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements,
	typename std::enable_if<
		ndim == 2 &&
		ETYPE == TransferElementType::SYMMETRIC_ISOTROPIC, int>::type*
)
{
	double initStart, initEnd;

	constexpr static double straightAngleRec = 1.0 / 180;

	initStart = eslog::time();
	eslog::info("       = LOOP TYPE                                                           MANUAL = \n");

	if (this->K == nullptr) {
		return loop<HeatTransferDataDescriptor, nodes, gps, ndim, edim, ETYPE>(action, ops, elements);
	}
	if ((action == ActionOperator::FILL) || (action == ActionOperator::SOLUTION))
	{
		return conditionsloop<DataDescriptor, nodes, gps, ndim, edim, ETYPE>(action, ops, interval, elements);
	}
	if (elements == 0) return measurements();

	typename DataDescriptor<nodes, gps, ndim, edim, ETYPE>::Element element;
	std::vector<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*> active; active.reserve(ops.size());

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if (elements > SIMD::size) {
				if ((*op)->isconst) {
					dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op)->simd(element);
				} else {
					active.push_back(dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op));
					active.back()->simd(element);
				}
			} else {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op)->peel(element, elements);
			}
		}
	}

	auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;

	const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
	bool rotateConductivity = mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	bool isIsotropic = mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	// it is dirty hack just to be sure that compiler must assume both variants (currently settings.sigma = 0 and diffusion_split = false)
	bool constConductivity = !settings.diffusion_split;
	bool constRotation = settings.sigma == 0;
	if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
		if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			if (ndim == 2) {
				rotateConductivity &= mat->coordinate_system.rotation.z.isset;
			}
			if (ndim == 3) {
				rotateConductivity &= mat->coordinate_system.rotation.x.isset | mat->coordinate_system.rotation.y.isset | mat->coordinate_system.rotation.z.isset;
			}
		}
	}

	bool hasHeatSource = false;
	bool constHeatSource = true;
	auto heatSourceEval = configuration.heat_source.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (heatSourceEval != configuration.heat_source.end()) {
		hasHeatSource = true;
		constHeatSource = heatSourceEval->second.evaluator != nullptr;
	}

	bool hasAdvection = false;
	bool constAdvection = true;
	auto advectionEval = configuration.translation_motions.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (advectionEval != configuration.translation_motions.end()) {
		hasAdvection = true;
		constAdvection = advectionEval->second.x.evaluator != nullptr && advectionEval->second.y.evaluator != nullptr && advectionEval->second.z.evaluator != nullptr;
	}

	bool cooToGP = mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN;
	bool isSpherical = mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::SPHERICAL;
	bool computeK = action == ActionOperator::ASSEMBLE || action == ActionOperator::REASSEMBLE;
	bool computeGradient = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.gradient;
	bool computeFlux = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.flux;
	bool computeConductivity = computeK | computeFlux;
	bool getTemp = computeGradient || computeFlux;
	bool isfullMatrix = this->K->shape == Matrix_Shape::FULL;

	initEnd = eslog::time();

	OutputParameterIterator stiffness(this->elements.stiffness, interval);

	if((action == ActionOperator::VOID) || (action == ActionOperator::SOLUTION))
	{
		std::cout<<"UNSUPPOERTED ACTION"<<std::endl;
		return measurements();
	}

	if( hasAdvection ||
		hasHeatSource ||
		getTemp)
	{
		std::cout<<"UNSUPPORTED OPERATOR"<<std::endl;
		return measurements();
	}

	if (!constRotation || !constConductivity)
	{
		std::cout<<"NON-CONST ROTATION OR CONDUCTIVITY"<<std::endl;
		return measurements();
	}


	// std::cout<<"ELEMENT TYPE IS "<<ETYPE<<std::endl;

	stiffness += SIMD::size;
	for (size_t s = 0; s < SIMD::size; ++s) {
		++procNodes;
	}

	// std::cout<<"WATCH OUT FOR THIS"<<std::endl;

	double start, end;
	if((action == ActionOperator::ASSEMBLE) || (action == ActionOperator::REASSEMBLE))
	{
		if(!cooToGP && !rotateConductivity && isIsotropic)
		{
			start = eslog::time();
			// std::cout<<"CARTESIAN ISOTROPIC NO ROTATION -- ISOTROPIC2D"<<std::endl;
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xFACE);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xCAFE);
			}
			esint chunks = elements / SIMD::size;
			for (esint c = 1; c < chunks; ++c) {

				for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
					for (size_t n = 0; n < nodes; ++n) {
						for (size_t d = 0; d < ndim; ++d) {
							element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
						}
					}
				}

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
				}

				double * __restrict__ out = stiffness.data;
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]) * element.conductivity[gp];
					for (size_t n = 0, i = 0; n < nodes; ++n) {
						SIMD nx = element.dND[gp][n][0];
						SIMD ny = element.dND[gp][n][1];
						for (size_t m = n; m < nodes; ++m, ++i) {
							SIMD mx = element.dND[gp][m][0];
							SIMD my = element.dND[gp][m][1];

							SIMD res = load(out + i * SIMD::size);
							res = res + scale * (nx * mx + ny * my);
							store(out + i * SIMD::size, res);
						}
					}
				}
				stiffness+= SIMD::size;

			}
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xDEAD);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xFADE);
			}
			end = eslog::time();
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
	return measurements(initEnd - initStart, end - start);
}

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t ETYPE>
Assembler::measurements HeatTransfer::manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements,
	typename std::enable_if<
		ndim == 2 &&
		ETYPE == TransferElementType::SYMMETRIC_GENERAL, int>::type*
)
{
	double initStart, initEnd;

	constexpr static double straightAngleRec = 1.0 / 180;

	initStart = eslog::time();
	eslog::info("       = LOOP TYPE                                                           MANUAL = \n");

	if (this->K == nullptr) {
		// std::cout<<"K IS NULL"<<std::endl;
		return loop<HeatTransferDataDescriptor, nodes, gps, ndim, edim, ETYPE>(action, ops, elements);
	}
	if ((action == ActionOperator::FILL) || (action == ActionOperator::SOLUTION))
	{
		return conditionsloop<DataDescriptor, nodes, gps, ndim, edim, ETYPE>(action, ops, interval, elements);
	}
	if (elements == 0) return measurements();

	typename DataDescriptor<nodes, gps, ndim, edim, ETYPE>::Element element;
	std::vector<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*> active; active.reserve(ops.size());

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if (elements > SIMD::size) {
				if ((*op)->isconst) {
					dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op)->simd(element);
				} else {
					active.push_back(dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op));
					active.back()->simd(element);
				}
			} else {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op)->peel(element, elements);
			}
		}
	}

	auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;

	const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
	bool rotateConductivity = mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	bool isIsotropic = mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	// it is dirty hack just to be sure that compiler must assume both variants (currently settings.sigma = 0 and diffusion_split = false)
	bool constConductivity = !settings.diffusion_split;
	bool constRotation = settings.sigma == 0;
	if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
		if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			if (ndim == 2) {
				rotateConductivity &= mat->coordinate_system.rotation.z.isset;
			}
			if (ndim == 3) {
				rotateConductivity &= mat->coordinate_system.rotation.x.isset | mat->coordinate_system.rotation.y.isset | mat->coordinate_system.rotation.z.isset;
			}
		}
	}

	bool hasHeatSource = false;
	bool constHeatSource = true;
	auto heatSourceEval = configuration.heat_source.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (heatSourceEval != configuration.heat_source.end()) {
		hasHeatSource = true;
		constHeatSource = heatSourceEval->second.evaluator != nullptr;
	}

	bool hasAdvection = false;
	bool constAdvection = true;
	auto advectionEval = configuration.translation_motions.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (advectionEval != configuration.translation_motions.end()) {
		hasAdvection = true;
		constAdvection = advectionEval->second.x.evaluator != nullptr && advectionEval->second.y.evaluator != nullptr && advectionEval->second.z.evaluator != nullptr;
	}

	bool cooToGP = mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN;
	bool isSpherical = mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::SPHERICAL;
	bool computeK = action == ActionOperator::ASSEMBLE || action == ActionOperator::REASSEMBLE;
	bool computeGradient = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.gradient;
	bool computeFlux = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.flux;
	bool computeConductivity = computeK | computeFlux;
	bool getTemp = computeGradient || computeFlux;
	bool isfullMatrix = this->K->shape == Matrix_Shape::FULL;

	initEnd = eslog::time();

	OutputParameterIterator stiffness(this->elements.stiffness, interval);

	if((action == ActionOperator::VOID) || (action == ActionOperator::SOLUTION) || (action == ActionOperator::FILL))
	{
		std::cout<<"UNSUPPOERTED ACTION"<<std::endl;
		return measurements();
	}

	if( hasAdvection ||
		hasHeatSource ||
		getTemp)
	{
		std::cout<<"UNSUPPORTED OPERATOR"<<std::endl;
		return measurements();
	}

	if (!constRotation || !constConductivity)
	{
		std::cout<<"NON-CONST ROTATION OR CONDUCTIVITY"<<std::endl;
		return measurements();
	}


	// std::cout<<"ELEMENT TYPE IS "<<ETYPE<<std::endl;

	stiffness += SIMD::size;
	for (size_t s = 0; s < SIMD::size; ++s) {
		++procNodes;
	}

	// std::cout<<"WATCH OUT FOR THIS"<<std::endl;

	double start, end;
	if((action == ActionOperator::ASSEMBLE) || (action == ActionOperator::REASSEMBLE))
	{
		if(!cooToGP && !rotateConductivity)
		{
			start = eslog::time();
			// std::cout<<"CARTESIAN NON-ISOTROPIC NO ROTATION -- DIAGONAL2D"<<std::endl;
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xFACE);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xCAFE);
			}
			esint chunks = elements / SIMD::size;
			for (esint c = 1; c < chunks; ++c) {
				// coo.simd(element);
				for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
					for (size_t n = 0; n < nodes; ++n) {
						for (size_t d = 0; d < ndim; ++d) {
							element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
						}
					}
				}

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
				}

				double * __restrict__ out = stiffness.data;
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD c00 = element.conductivity[gp][0];
					SIMD c10 = element.conductivity[gp][1], c11 = element.conductivity[gp][2];
					SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]);
					for (size_t n = 0, i = 0; n < nodes; ++n) {
						SIMD nx = element.dND[gp][n][0];
						SIMD ny = element.dND[gp][n][1];
						SIMD a = nx * c00 + ny * c10;
						SIMD b = nx * c10 + ny * c11;
						for (size_t m = n; m < nodes; ++m, ++i) {
							SIMD mx = element.dND[gp][m][0];
							SIMD my = element.dND[gp][m][1];

							SIMD res = load(out + i * SIMD::size);
							res = res + scale * (a * mx + b * my);
							store(out + i * SIMD::size, res);
						}
					}
				}
				stiffness+= SIMD::size;
			}
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xDEAD);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xFADE);
			}
			end = eslog::time();
		}
		if(!cooToGP && rotateConductivity)
		{
			start = eslog::time();
			// std::cout<<"CARTESIAN NON-ISOTROPIC YES ROTATION SYMETRIC -- SYMMETRIC2D,CARTESIAN2D"<<std::endl;
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xFACE);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xCAFE);
			}
			esint chunks = elements / SIMD::size;
			for (esint c = 1; c < chunks; ++c) {
				for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
					for (size_t n = 0; n < nodes; ++n) {
						for (size_t d = 0; d < ndim; ++d) {
							element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
						}
					}
				}

				//integration.simd(element);
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
				}

				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD angle = element.ecf.center[gp][0];
					for (size_t s = 0; s < SIMD::size; ++s) {
						element.cossin[gp][0][s] = std::cos(M_PI * angle[s] * straightAngleRec);
						element.cossin[gp][1][s] = std::sin(M_PI * angle[s] * straightAngleRec);
					}
				}

				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD c00 = element.ecf.conductivity[gp][0];
					SIMD c01 = element.ecf.conductivity[gp][1], c11 = element.ecf.conductivity[gp][2];
					SIMD cos = element.cossin[gp][0];
					SIMD sin = element.cossin[gp][1];

					element.conductivity[gp][0] = (cos * c00 - sin * c01) * cos - (cos * c01 - sin * c11) * sin;
					element.conductivity[gp][1] = (cos * c00 - sin * c01) * sin + (cos * c01 - sin * c11) * cos;
					element.conductivity[gp][2] = (sin * c00 + cos * c01) * sin + (sin * c01 + cos * c11) * cos;
				}

				double * __restrict__ out = stiffness.data;
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD c00 = element.conductivity[gp][0];
					SIMD c10 = element.conductivity[gp][1], c11 = element.conductivity[gp][2];
					SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]);
					for (size_t n = 0, i = 0; n < nodes; ++n) {
						SIMD nx = element.dND[gp][n][0];
						SIMD ny = element.dND[gp][n][1];
						SIMD a = nx * c00 + ny * c10;
						SIMD b = nx * c10 + ny * c11;
						for (size_t m = n; m < nodes; ++m, ++i) {
							SIMD mx = element.dND[gp][m][0];
							SIMD my = element.dND[gp][m][1];

							SIMD res = load(out + i * SIMD::size);
							res = res + scale * (a * mx + b * my);
							store(out + i * SIMD::size, res);
						}
					}
				}
				stiffness+= SIMD::size;

			}
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xDEAD);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xFADE);
			}
			end = eslog::time();
		}
		if(cooToGP && rotateConductivity)
			{
				start = eslog::time();
				// std::cout<<"NON-CARTESIAN NON-ISOTROPIC YES ROTATION SYMMETRIC -- CYLINDRICAL2D"<<std::endl;
				if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xFACE);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xCAFE);
			}
				esint chunks = elements / SIMD::size;
				for (esint c = 1; c < chunks; ++c) {
					// cooToGP.simd(element);
					for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
						for (size_t n = 0; n < nodes; ++n) {
							for (size_t d = 0; d < ndim; ++d) {
								element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
							}
						}
					}
					for (size_t gp = 0; gp < gps; ++gp) {
						for (size_t d = 0; d < ndim; ++d) {
							element.gpcoords[gp][d] = zeros();
						}
						for (size_t n = 0; n < nodes; ++n) {
							for (size_t d = 0; d < ndim; ++d) {
								element.gpcoords[gp][d] = element.gpcoords[gp][d] + load1(element.N[gp][n]) * element.coords[n][d];
							}
						}
					}

					//integration.simd(element);
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
					}

					for (size_t gp = 0; gp < gps; ++gp) {
						SIMD cooX =    element.gpcoords[gp][0];
						SIMD cooY =    element.gpcoords[gp][1];
						SIMD centerX = element.ecf.center[gp][0];
						SIMD centerY = element.ecf.center[gp][1];
						SIMD distanceX = cooX - centerX;
						SIMD distanceY = cooY - centerY;
						for (size_t s = 0; s < SIMD::size; ++s) {
							double rot = std::atan2(distanceY[s], distanceX[s]);
							element.cossin[gp][0][s] = std::cos(rot);
							element.cossin[gp][1][s] = std::sin(rot);
						}
					}

					for (size_t gp = 0; gp < gps; ++gp) {
						SIMD c00 = element.ecf.conductivity[gp][0];
						SIMD c01 = element.ecf.conductivity[gp][1], c11 = element.ecf.conductivity[gp][2];
						SIMD cos = element.cossin[gp][0];
						SIMD sin = element.cossin[gp][1];

						element.conductivity[gp][0] = (cos * c00 - sin * c01) * cos - (cos * c01 - sin * c11) * sin;
						element.conductivity[gp][1] = (cos * c00 - sin * c01) * sin + (cos * c01 - sin * c11) * cos;
						element.conductivity[gp][2] = (sin * c00 + cos * c01) * sin + (sin * c01 + cos * c11) * cos;
					}

					double * __restrict__ out = stiffness.data;
					for (size_t gp = 0; gp < gps; ++gp) {
						SIMD c00 = element.conductivity[gp][0];
						SIMD c10 = element.conductivity[gp][1], c11 = element.conductivity[gp][2];
						SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]);
						for (size_t n = 0, i = 0; n < nodes; ++n) {
							SIMD nx = element.dND[gp][n][0];
							SIMD ny = element.dND[gp][n][1];
							SIMD a = nx * c00 + ny * c10;
							SIMD b = nx * c10 + ny * c11;
							for (size_t m = n; m < nodes; ++m, ++i) {
								SIMD mx = element.dND[gp][m][0];
								SIMD my = element.dND[gp][m][1];

								SIMD res = load(out + i * SIMD::size);
								res = res + scale * (a * mx + b * my);
								store(out + i * SIMD::size, res);
							}
						}
					}
					stiffness+= SIMD::size;
				}
				if (action == ActionOperator::ASSEMBLE)
				{
					__SSC_MARK(0xDEAD);
				}
				if (action == ActionOperator::REASSEMBLE)
				{
					__SSC_MARK(0xFADE);
				}
				end = eslog::time();
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
	return measurements(initEnd - initStart, end - start);
}

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t ETYPE>
Assembler::measurements HeatTransfer::manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements,
	typename std::enable_if<
		ndim == 2 &&
		ETYPE == TransferElementType::ASYMMETRIC_ISOTROPIC, int>::type*
)
{
	return measurements();
}

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t ETYPE>
Assembler::measurements HeatTransfer::manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements,
	typename std::enable_if<
		ndim == 2 &&
		ETYPE == TransferElementType::ASYMMETRIC_GENERAL, int>::type*
)
{
	double initStart, initEnd;

	constexpr static double straightAngleRec = 1.0 / 180;

	initStart = eslog::time();
	eslog::info("       = LOOP TYPE                                                           MANUAL = \n");

	if (this->K == nullptr) {
		std::cout<<"K IS NULL"<<std::endl;
		return loop<HeatTransferDataDescriptor, nodes, gps, ndim, edim, ETYPE>(action, ops, elements);
	}
	if ((action == ActionOperator::FILL) || (action == ActionOperator::SOLUTION))
	{
		return conditionsloop<DataDescriptor, nodes, gps, ndim, edim, ETYPE>(action, ops, interval, elements);
	}
	if (elements == 0) return measurements();

	typename DataDescriptor<nodes, gps, ndim, edim, ETYPE>::Element element;
	std::vector<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*> active; active.reserve(ops.size());

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if (elements > SIMD::size) {
				if ((*op)->isconst) {
					dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op)->simd(element);
				} else {
					active.push_back(dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op));
					active.back()->simd(element);
				}
			} else {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op)->peel(element, elements);
			}
		}
	}

	auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;

	const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
	bool rotateConductivity = mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	bool isIsotropic = mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	// it is dirty hack just to be sure that compiler must assume both variants (currently settings.sigma = 0 and diffusion_split = false)
	bool constConductivity = !settings.diffusion_split;
	bool constRotation = settings.sigma == 0;
	if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
		if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			if (ndim == 2) {
				rotateConductivity &= mat->coordinate_system.rotation.z.isset;
			}
			if (ndim == 3) {
				rotateConductivity &= mat->coordinate_system.rotation.x.isset | mat->coordinate_system.rotation.y.isset | mat->coordinate_system.rotation.z.isset;
			}
		}
	}

	bool hasHeatSource = false;
	bool constHeatSource = true;
	auto heatSourceEval = configuration.heat_source.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (heatSourceEval != configuration.heat_source.end()) {
		hasHeatSource = true;
		constHeatSource = heatSourceEval->second.evaluator != nullptr;
	}

	bool hasAdvection = false;
	bool constAdvection = true;
	auto advectionEval = configuration.translation_motions.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (advectionEval != configuration.translation_motions.end()) {
		hasAdvection = true;
		constAdvection = advectionEval->second.x.evaluator != nullptr && advectionEval->second.y.evaluator != nullptr && advectionEval->second.z.evaluator != nullptr;
	}

	bool cooToGP = mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN;
	bool isSpherical = mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::SPHERICAL;
	bool computeK = action == ActionOperator::ASSEMBLE || action == ActionOperator::REASSEMBLE;
	bool computeGradient = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.gradient;
	bool computeFlux = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.flux;
	bool computeConductivity = computeK | computeFlux;
	bool getTemp = computeGradient || computeFlux;
	bool isfullMatrix = this->K->shape == Matrix_Shape::FULL;

	initEnd = eslog::time();

	OutputParameterIterator stiffness(this->elements.stiffness, interval);

	if((action == ActionOperator::VOID) || (action == ActionOperator::SOLUTION) || (action == ActionOperator::FILL))
	{
		std::cout<<"UNSUPPOERTED ACTION"<<std::endl;
		return measurements();
	}

	if( hasAdvection ||
		hasHeatSource ||
		getTemp)
	{
		std::cout<<"UNSUPPORTED OPERATOR"<<std::endl;
		return measurements();
	}

	if (!constRotation || !constConductivity)
	{
		std::cout<<"NON-CONST ROTATION OR CONDUCTIVITY"<<std::endl;
		return measurements();
	}


	// std::cout<<"ELEMENT TYPE IS "<<ETYPE<<std::endl;

	stiffness += SIMD::size;
	for (size_t s = 0; s < SIMD::size; ++s) {
		++procNodes;
	}

	// std::cout<<"WATCH OUT FOR THIS"<<std::endl;

	double start, end;
	if((action == ActionOperator::ASSEMBLE) || (action == ActionOperator::REASSEMBLE))
	{
		if(!cooToGP && rotateConductivity)
		{
			start = eslog::time();
			// std::cout<<"CARTESIAN NON-ISOTROPIC YES ROTATION NON-SYMETRIC -- ANISOTROPIC2D"<<std::endl;
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xFACE);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xCAFE);
			}
			esint chunks = elements / SIMD::size;
			for (esint c = 1; c < chunks; ++c) {
				for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
					for (size_t n = 0; n < nodes; ++n) {
						for (size_t d = 0; d < ndim; ++d) {
							element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
						}
					}
				}

				//integration.simd(element);
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
				}

				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD angle = element.ecf.center[gp][0];
					for (size_t s = 0; s < SIMD::size; ++s) {
						element.cossin[gp][0][s] = std::cos(M_PI * angle[s] * straightAngleRec);
						element.cossin[gp][1][s] = std::sin(M_PI * angle[s] * straightAngleRec);
					}
				}

				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD origin0 = element.ecf.conductivity[gp][0];
					SIMD origin1 = element.ecf.conductivity[gp][1];
					SIMD origin2 = element.ecf.conductivity[gp][2];
					SIMD origin3 = element.ecf.conductivity[gp][3];
					SIMD cos = element.cossin[gp][0];
					SIMD sin = element.cossin[gp][1];

					element.conductivity[gp][0] = (cos * origin0 - sin * origin2) * cos - (cos * origin1 - sin * origin3) * sin;
					element.conductivity[gp][1] = (cos * origin0 - sin * origin2) * sin + (cos * origin1 - sin * origin3) * cos;
					element.conductivity[gp][2] = (sin * origin0 + cos * origin2) * cos - (sin * origin1 + cos * origin3) * sin;
					element.conductivity[gp][3] = (sin * origin0 + cos * origin2) * sin + (sin * origin1 + cos * origin3) * cos;
				}

				double * __restrict__ out = stiffness.data;
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD c00 = element.conductivity[gp][0], c01 = element.conductivity[gp][2];
					SIMD c10 = element.conductivity[gp][1], c11 = element.conductivity[gp][3];
					SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]);
					for (size_t n = 0; n < nodes; ++n) {
						SIMD nx = element.dND[gp][n][0];
						SIMD ny = element.dND[gp][n][1];
						SIMD a = nx * c00 + ny * c01;
						SIMD b = nx * c10 + ny * c11;
						for (size_t m = 0; m < nodes; ++m) {
							SIMD mx = element.dND[gp][m][0];
							SIMD my = element.dND[gp][m][1];
							SIMD res = load(out + (n * nodes + m) * SIMD::size);
							res = res + scale * (a * mx + b * my);
							store(out + (n * nodes + m) * SIMD::size, res);
						}
					}
				}
				stiffness+= SIMD::size;
			}
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xDEAD);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xFADE);
			}
			end = eslog::time();
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
	return measurements(initEnd - initStart, end - start);
}

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t ETYPE>
Assembler::measurements HeatTransfer::manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements,
	typename std::enable_if<
		ndim == 3 &&
		ETYPE == TransferElementType::SYMMETRIC_ISOTROPIC, int>::type*
)
{
	double initStart, initEnd;

	constexpr static double straightAngleRec = 1.0 / 180;

	initStart = eslog::time();
	eslog::info("       = LOOP TYPE                                                           MANUAL = \n");

	if (this->K == nullptr) {
		std::cout<<"K IS NULL"<<std::endl;
		return loop<HeatTransferDataDescriptor, nodes, gps, ndim, edim, ETYPE>(action, ops, elements);
	}
	if ((action == ActionOperator::FILL) || (action == ActionOperator::SOLUTION))
	{
		return conditionsloop<DataDescriptor, nodes, gps, ndim, edim, ETYPE>(action, ops, interval, elements);
	}
	if (elements == 0) return measurements();

	typename DataDescriptor<nodes, gps, ndim, edim, ETYPE>::Element element;
	std::vector<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*> active; active.reserve(ops.size());

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if (elements > SIMD::size) {
				if ((*op)->isconst) {
					dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op)->simd(element);
				} else {
					active.push_back(dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op));
					active.back()->simd(element);
				}
			} else {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op)->peel(element, elements);
			}
		}
	}

	auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;

	const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
	bool rotateConductivity = mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	bool isIsotropic = mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	// it is dirty hack just to be sure that compiler must assume both variants (currently settings.sigma = 0 and diffusion_split = false)
	bool constConductivity = !settings.diffusion_split;
	bool constRotation = settings.sigma == 0;
	if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
		if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			if (ndim == 2) {
				rotateConductivity &= mat->coordinate_system.rotation.z.isset;
			}
			if (ndim == 3) {
				rotateConductivity &= mat->coordinate_system.rotation.x.isset | mat->coordinate_system.rotation.y.isset | mat->coordinate_system.rotation.z.isset;
			}
		}
	}

	bool hasHeatSource = false;
	bool constHeatSource = true;
	auto heatSourceEval = configuration.heat_source.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (heatSourceEval != configuration.heat_source.end()) {
		hasHeatSource = true;
		constHeatSource = heatSourceEval->second.evaluator != nullptr;
	}

	bool hasAdvection = false;
	bool constAdvection = true;
	auto advectionEval = configuration.translation_motions.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (advectionEval != configuration.translation_motions.end()) {
		hasAdvection = true;
		constAdvection = advectionEval->second.x.evaluator != nullptr && advectionEval->second.y.evaluator != nullptr && advectionEval->second.z.evaluator != nullptr;
	}

	bool cooToGP = mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN;
	bool isSpherical = mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::SPHERICAL;
	bool computeK = action == ActionOperator::ASSEMBLE || action == ActionOperator::REASSEMBLE;
	bool computeGradient = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.gradient;
	bool computeFlux = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.flux;
	bool computeConductivity = computeK | computeFlux;
	bool getTemp = computeGradient || computeFlux;
	bool isfullMatrix = this->K->shape == Matrix_Shape::FULL;

	initEnd = eslog::time();

	OutputParameterIterator stiffness(this->elements.stiffness, interval);

	if((action == ActionOperator::VOID) || (action == ActionOperator::SOLUTION) || (action == ActionOperator::FILL))
	{
		std::cout<<"UNSUPPOERTED ACTION"<<std::endl;
		return measurements();
	}

	if( hasAdvection ||
		hasHeatSource ||
		getTemp)
	{
		std::cout<<"UNSUPPORTED OPERATOR"<<std::endl;
		return measurements();
	}

	if (!constRotation || !constConductivity)
	{
		std::cout<<"NON-CONST ROTATION OR CONDUCTIVITY"<<std::endl;
		return measurements();
	}

	// std::cout<<"ELEMENT TYPE IS "<<ETYPE<<std::endl;

	stiffness += SIMD::size;
	for (size_t s = 0; s < SIMD::size; ++s) {
		++procNodes;
	}

	// std::cout<<"WATCH OUT FOR THIS"<<std::endl;

	double start, end;
	if((action == ActionOperator::ASSEMBLE) || (action == ActionOperator::REASSEMBLE))
	{
		if(!cooToGP && !rotateConductivity && isIsotropic)
		{
			start = eslog::time();
			// std::cout<<"CARTESIAN ISOTROPIC NO ROTATION -- ISOTROPIC3D"<<std::endl;
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xFACE);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xCAFE);
			}
			esint chunks = elements / SIMD::size;
			for (esint c = 1; c < chunks; ++c) {
				//coo
				for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
					for (size_t n = 0; n < nodes; ++n) {
						for (size_t d = 0; d < ndim; ++d) {
							element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
						}
					}
				}
				//integrate
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD jacobian0 = zeros(), jacobian1 = zeros(), jacobian2 = zeros(), jacobian3 = zeros(), jacobian4 = zeros();
					SIMD jacobian5 = zeros(), jacobian6 = zeros(), jacobian7 = zeros(), jacobian8 = zeros();

					for (size_t n = 0; n < nodes; ++n) {
						SIMD coordsX = element.coords[n][0];
						SIMD coordsY = element.coords[n][1];
						SIMD coordsZ = element.coords[n][2];
						SIMD dNX = load1(element.dN[gp][n][0]);
						SIMD dNY = load1(element.dN[gp][n][1]);
						SIMD dNZ = load1(element.dN[gp][n][2]);

						jacobian0 = jacobian0 + dNX * coordsX;
						jacobian1 = jacobian1 + dNX * coordsY;
						jacobian2 = jacobian2 + dNX * coordsZ;
						jacobian3 = jacobian3 + dNY * coordsX;
						jacobian4 = jacobian4 + dNY * coordsY;
						jacobian5 = jacobian5 + dNY * coordsZ;
						jacobian6 = jacobian6 + dNZ * coordsX;
						jacobian7 = jacobian7 + dNZ * coordsY;
						jacobian8 = jacobian8 + dNZ * coordsZ;
					}

					element.det[gp] =
							+ jacobian0 * jacobian4 * jacobian8
							+ jacobian1 * jacobian5 * jacobian6
							+ jacobian2 * jacobian3 * jacobian7
							- jacobian2 * jacobian4 * jacobian6
							- jacobian1 * jacobian3 * jacobian8
							- jacobian0 * jacobian5 * jacobian7;

					SIMD detJx = ones() / element.det[gp];
					SIMD inv0 = detJx * ( jacobian8 * jacobian4 - jacobian7 * jacobian5);
					SIMD inv1 = detJx * (-jacobian8 * jacobian1 + jacobian7 * jacobian2);
					SIMD inv2 = detJx * ( jacobian5 * jacobian1 - jacobian4 * jacobian2);
					SIMD inv3 = detJx * (-jacobian8 * jacobian3 + jacobian6 * jacobian5);
					SIMD inv4 = detJx * ( jacobian8 * jacobian0 - jacobian6 * jacobian2);
					SIMD inv5 = detJx * (-jacobian5 * jacobian0 + jacobian3 * jacobian2);
					SIMD inv6 = detJx * ( jacobian7 * jacobian3 - jacobian6 * jacobian4);
					SIMD inv7 = detJx * (-jacobian7 * jacobian0 + jacobian6 * jacobian1);
					SIMD inv8 = detJx * ( jacobian4 * jacobian0 - jacobian3 * jacobian1);

					for (size_t n = 0; n < nodes; ++n) {
						SIMD dNX = load1(element.dN[gp][n][0]);
						SIMD dNY = load1(element.dN[gp][n][1]);
						SIMD dNZ = load1(element.dN[gp][n][2]);
						element.dND[gp][n][0] = inv0 * dNX + inv1 * dNY + inv2 * dNZ;
						element.dND[gp][n][1] = inv3 * dNX + inv4 * dNY + inv5 * dNZ;
						element.dND[gp][n][2] = inv6 * dNX + inv7 * dNY + inv8 * dNZ;
					}
				}
				//stiffness symmetric isotropic
				double * __restrict__ out = stiffness.data;
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD scale = element.det[gp] * load1(element.w[gp]) * element.conductivity[gp];
					for (size_t n = 0, i = 0; n < nodes; ++n) {
						SIMD nx = element.dND[gp][n][0];
						SIMD ny = element.dND[gp][n][1];
						SIMD nz = element.dND[gp][n][2];
						for (size_t m = n; m < nodes; ++m, ++i) {
							SIMD mx = element.dND[gp][m][0];
							SIMD my = element.dND[gp][m][1];
							SIMD mz = element.dND[gp][m][2];

							SIMD res = load(out + i * SIMD::size);
							res = res + scale * (nx * mx + ny * my + nz * mz);
							store(out + i * SIMD::size, res);
						}
					}
				}
				stiffness+= SIMD::size;
			}
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xDEAD);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xFADE);
			}
			end = eslog::time();
		}
	}


	// if(action == ActionOperator::REASSEMBLE)
	// {
	// 	start = eslog::time();
	// 	__SSC_MARK(0xCAFE);
	// 	esint chunks = elements / SIMD::size;
	// 	for (esint c = 1; c < chunks; ++c) {
	// 		if (cooToGP) {
	// 			cooAndGps.simd(element);
	// 		} else {
	// 			coo.simd(element);
	// 		}
	// 		integration.simd(element);
	// 		if (getTemp) {
	// 			temp.simd(element);
	// 		}
	// 		if (computeConductivity) {
	// 			if (!constConductivity) {
	// 				updateConductivity<DataDescriptor, nodes, gps, ndim, edim, ETYPE>()(element, mat);
	// 			}
	// 			if (rotateConductivity) {
	// 				if (!constRotation) {
	// 					updateRotation<DataDescriptor, nodes, gps, ndim, edim, ETYPE>()(element, mat);
	// 				}
	// 				rotation(element, mat);
	// 			}
	// 		}
	//
	// 		if (computeK) {
	// 			if (hasAdvection) {
	// 				if (!constAdvection) {
	// 					updateTM(element, advectionEval->second.x.evaluator);
	// 				}
	// 				advection.simd(element);
	// 			}
	// 			stiffness.simd(element);
	// 			if (hasHeatSource) {
	// 				if (!constHeatSource) {
	// 					updateHS(element, heatSourceEval->second.evaluator);
	// 				}
	// 				heatSource.simd(element);
	// 			}
	// 		}
	// 		if (action == ActionOperator::FILL) {
	// 			if (isfullMatrix) {
	// 				fullFiller.simd(element);
	// 			} else {
	// 				upperFiller.simd(element);
	// 			}
	// 			rhsFiller.simd(element);
	// 		}
	// 		if (computeGradient) {
	// 			gradient.simd(element);
	// 		}
	// 		if (computeFlux) {
	// 			flux.simd(element);
	// 		}
	// 	}
	// 	__SSC_MARK(0xFADE);
	// 	end = eslog::time();
	// }


	// if(action == ActionOperator::REASSEMBLE)
	// {
	// 	start = eslog::time();
	// 	__SSC_MARK(0xCAFE);
	// 	esint chunks = elements / SIMD::size;
	// 	for (esint c = 1; c < chunks; ++c) {
	// 		if (cooToGP) {
	// 			cooAndGps.simd(element);
	// 		} else {
	// 			coo.simd(element);
	// 		}
	// 		integration.simd(element);
	// 		if (getTemp) {
	// 			temp.simd(element);
	// 		}
	// 		if (computeConductivity) {
	// 			if (!constConductivity) {
	// 				updateConductivity<DataDescriptor, nodes, gps, ndim, edim, ETYPE>()(element, mat);
	// 			}
	// 			if (rotateConductivity) {
	// 				if (!constRotation) {
	// 					updateRotation<DataDescriptor, nodes, gps, ndim, edim, ETYPE>()(element, mat);
	// 				}
	// 				rotation(element, mat);
	// 			}
	// 		}
	//
	// 		if (computeK) {
	// 			if (hasAdvection) {
	// 				if (!constAdvection) {
	// 					updateTM(element, advectionEval->second.x.evaluator);
	// 				}
	// 				advection.simd(element);
	// 			}
	// 			stiffness.simd(element);
	// 			if (hasHeatSource) {
	// 				if (!constHeatSource) {
	// 					updateHS(element, heatSourceEval->second.evaluator);
	// 				}
	// 				heatSource.simd(element);
	// 			}
	// 		}
	// 		if (action == ActionOperator::FILL) {
	// 			if (isfullMatrix) {
	// 				fullFiller.simd(element);
	// 			} else {
	// 				upperFiller.simd(element);
	// 			}
	// 			rhsFiller.simd(element);
	// 		}
	// 		if (computeGradient) {
	// 			gradient.simd(element);
	// 		}
	// 		if (computeFlux) {
	// 			flux.simd(element);
	// 		}
	// 	}
	// 	__SSC_MARK(0xFADE);
	// 	end = eslog::time();
	// }


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
	return measurements(initEnd - initStart, end - start);
}

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t ETYPE>
Assembler::measurements HeatTransfer::manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements,
	typename std::enable_if<
		ndim == 3 &&
		ETYPE == TransferElementType::SYMMETRIC_GENERAL, int>::type*
)
{
	double initStart, initEnd;

	constexpr static double straightAngleRec = 1.0 / 180;

	initStart = eslog::time();
	eslog::info("       = LOOP TYPE                                                           MANUAL = \n");

	if (this->K == nullptr) {
		std::cout<<"K IS NULL"<<std::endl;
		return loop<HeatTransferDataDescriptor, nodes, gps, ndim, edim, ETYPE>(action, ops, elements);
	}
	if ((action == ActionOperator::FILL) || (action == ActionOperator::SOLUTION))
	{
		return conditionsloop<DataDescriptor, nodes, gps, ndim, edim, ETYPE>(action, ops, interval, elements);
	}
	if (elements == 0) return measurements();

	typename DataDescriptor<nodes, gps, ndim, edim, ETYPE>::Element element;
	std::vector<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*> active; active.reserve(ops.size());

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if (elements > SIMD::size) {
				if ((*op)->isconst) {
					dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op)->simd(element);
				} else {
					active.push_back(dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op));
					active.back()->simd(element);
				}
			} else {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op)->peel(element, elements);
			}
		}
	}

	auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;

	const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
	bool rotateConductivity = mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	bool isIsotropic = mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	// it is dirty hack just to be sure that compiler must assume both variants (currently settings.sigma = 0 and diffusion_split = false)
	bool constConductivity = !settings.diffusion_split;
	bool constRotation = settings.sigma == 0;
	if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
		if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			if (ndim == 2) {
				rotateConductivity &= mat->coordinate_system.rotation.z.isset;
			}
			if (ndim == 3) {
				rotateConductivity &= mat->coordinate_system.rotation.x.isset | mat->coordinate_system.rotation.y.isset | mat->coordinate_system.rotation.z.isset;
			}
		}
	}

	bool hasHeatSource = false;
	bool constHeatSource = true;
	auto heatSourceEval = configuration.heat_source.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (heatSourceEval != configuration.heat_source.end()) {
		hasHeatSource = true;
		constHeatSource = heatSourceEval->second.evaluator != nullptr;
	}

	bool hasAdvection = false;
	bool constAdvection = true;
	auto advectionEval = configuration.translation_motions.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (advectionEval != configuration.translation_motions.end()) {
		hasAdvection = true;
		constAdvection = advectionEval->second.x.evaluator != nullptr && advectionEval->second.y.evaluator != nullptr && advectionEval->second.z.evaluator != nullptr;
	}

	bool cooToGP = mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN;
	bool isSpherical = mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::SPHERICAL;
	bool computeK = action == ActionOperator::ASSEMBLE || action == ActionOperator::REASSEMBLE;
	bool computeGradient = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.gradient;
	bool computeFlux = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.flux;
	bool computeConductivity = computeK | computeFlux;
	bool getTemp = computeGradient || computeFlux;
	bool isfullMatrix = this->K->shape == Matrix_Shape::FULL;

	initEnd = eslog::time();

	OutputParameterIterator stiffness(this->elements.stiffness, interval);

	if((action == ActionOperator::VOID) || (action == ActionOperator::SOLUTION) || (action == ActionOperator::FILL))
	{
		std::cout<<"UNSUPPOERTED ACTION"<<std::endl;
		return measurements();
	}

	if( hasAdvection ||
		hasHeatSource ||
		getTemp)
	{
		std::cout<<"UNSUPPORTED OPERATOR"<<std::endl;
		return measurements();
	}

	if (!constRotation || !constConductivity)
	{
		std::cout<<"NON-CONST ROTATION OR CONDUCTIVITY"<<std::endl;
		return measurements();
	}

	// std::cout<<"ELEMENT TYPE IS "<<ETYPE<<std::endl;

	stiffness += SIMD::size;
	for (size_t s = 0; s < SIMD::size; ++s) {
		++procNodes;
	}

	// std::cout<<"WATCH OUT FOR THIS"<<std::endl;

	double start, end;
	if((action == ActionOperator::ASSEMBLE) || (action == ActionOperator::REASSEMBLE))
	{

		if(!cooToGP && !rotateConductivity && !isIsotropic)
		{
			start = eslog::time();
			// std::cout<<"CARTESIAN NON-ISOTROPIC NO ROTATION -- DIAGONAL3D"<<std::endl;
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xFACE);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xCAFE);
			}
			esint chunks = elements / SIMD::size;
			for (esint c = 1; c < chunks; ++c) {
				// coo.simd(element);
				for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
					for (size_t n = 0; n < nodes; ++n) {
						for (size_t d = 0; d < ndim; ++d) {
							element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
						}
					}
				}
				//integrate
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD jacobian0 = zeros(), jacobian1 = zeros(), jacobian2 = zeros(), jacobian3 = zeros(), jacobian4 = zeros();
					SIMD jacobian5 = zeros(), jacobian6 = zeros(), jacobian7 = zeros(), jacobian8 = zeros();

					for (size_t n = 0; n < nodes; ++n) {
						SIMD coordsX = element.coords[n][0];
						SIMD coordsY = element.coords[n][1];
						SIMD coordsZ = element.coords[n][2];
						SIMD dNX = load1(element.dN[gp][n][0]);
						SIMD dNY = load1(element.dN[gp][n][1]);
						SIMD dNZ = load1(element.dN[gp][n][2]);

						jacobian0 = jacobian0 + dNX * coordsX;
						jacobian1 = jacobian1 + dNX * coordsY;
						jacobian2 = jacobian2 + dNX * coordsZ;
						jacobian3 = jacobian3 + dNY * coordsX;
						jacobian4 = jacobian4 + dNY * coordsY;
						jacobian5 = jacobian5 + dNY * coordsZ;
						jacobian6 = jacobian6 + dNZ * coordsX;
						jacobian7 = jacobian7 + dNZ * coordsY;
						jacobian8 = jacobian8 + dNZ * coordsZ;
					}

					element.det[gp] =
							+ jacobian0 * jacobian4 * jacobian8
							+ jacobian1 * jacobian5 * jacobian6
							+ jacobian2 * jacobian3 * jacobian7
							- jacobian2 * jacobian4 * jacobian6
							- jacobian1 * jacobian3 * jacobian8
							- jacobian0 * jacobian5 * jacobian7;

					SIMD detJx = ones() / element.det[gp];
					SIMD inv0 = detJx * ( jacobian8 * jacobian4 - jacobian7 * jacobian5);
					SIMD inv1 = detJx * (-jacobian8 * jacobian1 + jacobian7 * jacobian2);
					SIMD inv2 = detJx * ( jacobian5 * jacobian1 - jacobian4 * jacobian2);
					SIMD inv3 = detJx * (-jacobian8 * jacobian3 + jacobian6 * jacobian5);
					SIMD inv4 = detJx * ( jacobian8 * jacobian0 - jacobian6 * jacobian2);
					SIMD inv5 = detJx * (-jacobian5 * jacobian0 + jacobian3 * jacobian2);
					SIMD inv6 = detJx * ( jacobian7 * jacobian3 - jacobian6 * jacobian4);
					SIMD inv7 = detJx * (-jacobian7 * jacobian0 + jacobian6 * jacobian1);
					SIMD inv8 = detJx * ( jacobian4 * jacobian0 - jacobian3 * jacobian1);

					for (size_t n = 0; n < nodes; ++n) {
						SIMD dNX = load1(element.dN[gp][n][0]);
						SIMD dNY = load1(element.dN[gp][n][1]);
						SIMD dNZ = load1(element.dN[gp][n][2]);
						element.dND[gp][n][0] = inv0 * dNX + inv1 * dNY + inv2 * dNZ;
						element.dND[gp][n][1] = inv3 * dNX + inv4 * dNY + inv5 * dNZ;
						element.dND[gp][n][2] = inv6 * dNX + inv7 * dNY + inv8 * dNZ;
					}
				}
				//stiffness symmetric non-isotropic
				double * __restrict__ out = stiffness.data;
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD scale = element.det[gp] * load1(element.w[gp]);
					SIMD c00 = element.conductivity[gp][0];
					SIMD c10 = element.conductivity[gp][1], c11 = element.conductivity[gp][3];
					SIMD c20 = element.conductivity[gp][2], c21 = element.conductivity[gp][4], c22 = element.conductivity[gp][5];
					for (size_t n = 0, i = 0; n < nodes; ++n) {
						SIMD nx = element.dND[gp][n][0];
						SIMD ny = element.dND[gp][n][1];
						SIMD nz = element.dND[gp][n][2];
						SIMD a = nx * c00 + ny * c10 + nz * c20;
						SIMD b = nx * c10 + ny * c11 + nz * c21;
						SIMD c = nx * c20 + ny * c21 + nz * c22;
						for (size_t m = n; m < nodes; ++m, ++i) {
							SIMD mx = element.dND[gp][m][0];
							SIMD my = element.dND[gp][m][1];
							SIMD mz = element.dND[gp][m][2];

							SIMD res = load(out + i * SIMD::size);
							res = res + scale * (a * mx + b * my + c * mz);
							store(out + i * SIMD::size, res);
						}
					}
				}
				stiffness+= SIMD::size;
			}
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xDEAD);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xFADE);
			}
			end = eslog::time();
		}
		if(!cooToGP && rotateConductivity)
		{
			start = eslog::time();
			// std::cout<<"CARTESIAN NON-ISOTROPIC YES ROTATION SYMETRIC -- SYMMETRIC3D,CARTESIAN3D"<<std::endl;
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xFACE);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xCAFE);
			}
			esint chunks = elements / SIMD::size;
			for (esint c = 1; c < chunks; ++c) {
				//coo
				for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
					for (size_t n = 0; n < nodes; ++n) {
						for (size_t d = 0; d < ndim; ++d) {
							element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
						}
					}
				}
				//integrate
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD jacobian0 = zeros(), jacobian1 = zeros(), jacobian2 = zeros(), jacobian3 = zeros(), jacobian4 = zeros();
					SIMD jacobian5 = zeros(), jacobian6 = zeros(), jacobian7 = zeros(), jacobian8 = zeros();

					for (size_t n = 0; n < nodes; ++n) {
						SIMD coordsX = element.coords[n][0];
						SIMD coordsY = element.coords[n][1];
						SIMD coordsZ = element.coords[n][2];
						SIMD dNX = load1(element.dN[gp][n][0]);
						SIMD dNY = load1(element.dN[gp][n][1]);
						SIMD dNZ = load1(element.dN[gp][n][2]);

						jacobian0 = jacobian0 + dNX * coordsX;
						jacobian1 = jacobian1 + dNX * coordsY;
						jacobian2 = jacobian2 + dNX * coordsZ;
						jacobian3 = jacobian3 + dNY * coordsX;
						jacobian4 = jacobian4 + dNY * coordsY;
						jacobian5 = jacobian5 + dNY * coordsZ;
						jacobian6 = jacobian6 + dNZ * coordsX;
						jacobian7 = jacobian7 + dNZ * coordsY;
						jacobian8 = jacobian8 + dNZ * coordsZ;
					}

					element.det[gp] =
							+ jacobian0 * jacobian4 * jacobian8
							+ jacobian1 * jacobian5 * jacobian6
							+ jacobian2 * jacobian3 * jacobian7
							- jacobian2 * jacobian4 * jacobian6
							- jacobian1 * jacobian3 * jacobian8
							- jacobian0 * jacobian5 * jacobian7;

					SIMD detJx = ones() / element.det[gp];
					SIMD inv0 = detJx * ( jacobian8 * jacobian4 - jacobian7 * jacobian5);
					SIMD inv1 = detJx * (-jacobian8 * jacobian1 + jacobian7 * jacobian2);
					SIMD inv2 = detJx * ( jacobian5 * jacobian1 - jacobian4 * jacobian2);
					SIMD inv3 = detJx * (-jacobian8 * jacobian3 + jacobian6 * jacobian5);
					SIMD inv4 = detJx * ( jacobian8 * jacobian0 - jacobian6 * jacobian2);
					SIMD inv5 = detJx * (-jacobian5 * jacobian0 + jacobian3 * jacobian2);
					SIMD inv6 = detJx * ( jacobian7 * jacobian3 - jacobian6 * jacobian4);
					SIMD inv7 = detJx * (-jacobian7 * jacobian0 + jacobian6 * jacobian1);
					SIMD inv8 = detJx * ( jacobian4 * jacobian0 - jacobian3 * jacobian1);

					for (size_t n = 0; n < nodes; ++n) {
						SIMD dNX = load1(element.dN[gp][n][0]);
						SIMD dNY = load1(element.dN[gp][n][1]);
						SIMD dNZ = load1(element.dN[gp][n][2]);
						element.dND[gp][n][0] = inv0 * dNX + inv1 * dNY + inv2 * dNZ;
						element.dND[gp][n][1] = inv3 * dNX + inv4 * dNY + inv5 * dNZ;
						element.dND[gp][n][2] = inv6 * dNX + inv7 * dNY + inv8 * dNZ;
					}
				}
				// cartesian
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD angleX = element.ecf.center[gp][0];
					SIMD angleY = element.ecf.center[gp][1];
					SIMD angleZ = element.ecf.center[gp][2];
					for (size_t s = 0; s < SIMD::size; ++s) {
						element.cossin[gp][0][s] = std::cos(M_PI * angleX[s] * straightAngleRec);
						element.cossin[gp][1][s] = std::cos(M_PI * angleY[s] * straightAngleRec);
						element.cossin[gp][2][s] = std::cos(M_PI * angleZ[s] * straightAngleRec);
						element.cossin[gp][3][s] = std::sin(M_PI * angleX[s] * straightAngleRec);
						element.cossin[gp][4][s] = std::sin(M_PI * angleY[s] * straightAngleRec);
						element.cossin[gp][5][s] = std::sin(M_PI * angleZ[s] * straightAngleRec);
					}
				}
				//rotate cartesian
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD cos0 = element.cossin[gp][0];
					SIMD cos1 = element.cossin[gp][1];
					SIMD cos2 = element.cossin[gp][2];
					SIMD sin0 = element.cossin[gp][3];
					SIMD sin1 = element.cossin[gp][4];
					SIMD sin2 = element.cossin[gp][5];

					SIMD t00 = cos1 * cos2;
					SIMD t01 = cos1 * sin2;
					SIMD t02 = -sin1;
					SIMD t10 = cos2 * sin0 * sin1 - cos0 * sin2;
					SIMD t11 = cos0 * cos2 + sin0 * sin1 * sin2;
					SIMD t12 = cos1 * sin0;
					SIMD t20 = sin0 * sin2 + cos0 * cos2 * sin1;
					SIMD t21 = cos0 * sin1 * sin2 - cos2 * sin0;
					SIMD t22 = cos0 * cos1;

					SIMD c00 = element.ecf.conductivity[gp][0];
					SIMD c01 = element.ecf.conductivity[gp][1], c11 = element.ecf.conductivity[gp][3];
					SIMD c02 = element.ecf.conductivity[gp][2], c12 = element.ecf.conductivity[gp][4], c22 = element.ecf.conductivity[gp][5];

					SIMD a = t00 * c00 + t10 * c01 + t20 * c02;
					SIMD b = t00 * c01 + t10 * c11 + t20 * c12;
					SIMD c = t00 * c02 + t10 * c12 + t20 * c22;
					element.conductivity[gp][0] = a * t00 + b * t10 + c * t20;
					element.conductivity[gp][1] = a * t01 + b * t11 + c * t21;
					element.conductivity[gp][2] = a * t02 + b * t12 + c * t22;

					a = t01 * c00 + t11 * c01 + t21 * c02;
					b = t01 * c01 + t11 * c11 + t21 * c12;
					c = t01 * c02 + t11 * c12 + t21 * c22;
					element.conductivity[gp][3] = a * t01 + b * t11 + c * t21;
					element.conductivity[gp][4] = a * t02 + b * t12 + c * t22;

					a = t02 * c00 + t12 * c01 + t22 * c02;
					b = t02 * c01 + t12 * c11 + t22 * c12;
					c = t02 * c02 + t12 * c12 + t22 * c22;
					element.conductivity[gp][5] = a * t02 + b * t12 + c * t22;
				}
				//stiffness symmetric non-isotropic
				double * __restrict__ out = stiffness.data;
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD scale = element.det[gp] * load1(element.w[gp]);
					SIMD c00 = element.conductivity[gp][0];
					SIMD c10 = element.conductivity[gp][1], c11 = element.conductivity[gp][3];
					SIMD c20 = element.conductivity[gp][2], c21 = element.conductivity[gp][4], c22 = element.conductivity[gp][5];
					for (size_t n = 0, i = 0; n < nodes; ++n) {
						SIMD nx = element.dND[gp][n][0];
						SIMD ny = element.dND[gp][n][1];
						SIMD nz = element.dND[gp][n][2];
						SIMD a = nx * c00 + ny * c10 + nz * c20;
						SIMD b = nx * c10 + ny * c11 + nz * c21;
						SIMD c = nx * c20 + ny * c21 + nz * c22;
						for (size_t m = n; m < nodes; ++m, ++i) {
							SIMD mx = element.dND[gp][m][0];
							SIMD my = element.dND[gp][m][1];
							SIMD mz = element.dND[gp][m][2];

							SIMD res = load(out + i * SIMD::size);
							res = res + scale * (a * mx + b * my + c * mz);
							store(out + i * SIMD::size, res);
						}
					}
				}
				stiffness+= SIMD::size;
			}
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xDEAD);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xFADE);
			}
			end = eslog::time();
		}
		if(cooToGP && rotateConductivity && !isSpherical)
		{
			start = eslog::time();
			// std::cout<<"NON-CARTESIAN NON-ISOTROPIC YES ROTATION SYMMETRIC NON-SPHERICAL -- CYLINDRICAL3D"<<std::endl;
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xFACE);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xCAFE);
			}
			esint chunks = elements / SIMD::size;
			for (esint c = 1; c < chunks; ++c) {
				// cooToGP.simd(element);
				for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
					for (size_t n = 0; n < nodes; ++n) {
						for (size_t d = 0; d < ndim; ++d) {
							element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
						}
					}
				}
				for (size_t gp = 0; gp < gps; ++gp) {
					for (size_t d = 0; d < ndim; ++d) {
						element.gpcoords[gp][d] = zeros();
					}
					for (size_t n = 0; n < nodes; ++n) {
						for (size_t d = 0; d < ndim; ++d) {
							element.gpcoords[gp][d] = element.gpcoords[gp][d] + load1(element.N[gp][n]) * element.coords[n][d];
						}
					}
				}
				//integrate
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD jacobian0 = zeros(), jacobian1 = zeros(), jacobian2 = zeros(), jacobian3 = zeros(), jacobian4 = zeros();
					SIMD jacobian5 = zeros(), jacobian6 = zeros(), jacobian7 = zeros(), jacobian8 = zeros();

					for (size_t n = 0; n < nodes; ++n) {
						SIMD coordsX = element.coords[n][0];
						SIMD coordsY = element.coords[n][1];
						SIMD coordsZ = element.coords[n][2];
						SIMD dNX = load1(element.dN[gp][n][0]);
						SIMD dNY = load1(element.dN[gp][n][1]);
						SIMD dNZ = load1(element.dN[gp][n][2]);

						jacobian0 = jacobian0 + dNX * coordsX;
						jacobian1 = jacobian1 + dNX * coordsY;
						jacobian2 = jacobian2 + dNX * coordsZ;
						jacobian3 = jacobian3 + dNY * coordsX;
						jacobian4 = jacobian4 + dNY * coordsY;
						jacobian5 = jacobian5 + dNY * coordsZ;
						jacobian6 = jacobian6 + dNZ * coordsX;
						jacobian7 = jacobian7 + dNZ * coordsY;
						jacobian8 = jacobian8 + dNZ * coordsZ;
					}

					element.det[gp] =
							+ jacobian0 * jacobian4 * jacobian8
							+ jacobian1 * jacobian5 * jacobian6
							+ jacobian2 * jacobian3 * jacobian7
							- jacobian2 * jacobian4 * jacobian6
							- jacobian1 * jacobian3 * jacobian8
							- jacobian0 * jacobian5 * jacobian7;

					SIMD detJx = ones() / element.det[gp];
					SIMD inv0 = detJx * ( jacobian8 * jacobian4 - jacobian7 * jacobian5);
					SIMD inv1 = detJx * (-jacobian8 * jacobian1 + jacobian7 * jacobian2);
					SIMD inv2 = detJx * ( jacobian5 * jacobian1 - jacobian4 * jacobian2);
					SIMD inv3 = detJx * (-jacobian8 * jacobian3 + jacobian6 * jacobian5);
					SIMD inv4 = detJx * ( jacobian8 * jacobian0 - jacobian6 * jacobian2);
					SIMD inv5 = detJx * (-jacobian5 * jacobian0 + jacobian3 * jacobian2);
					SIMD inv6 = detJx * ( jacobian7 * jacobian3 - jacobian6 * jacobian4);
					SIMD inv7 = detJx * (-jacobian7 * jacobian0 + jacobian6 * jacobian1);
					SIMD inv8 = detJx * ( jacobian4 * jacobian0 - jacobian3 * jacobian1);

					for (size_t n = 0; n < nodes; ++n) {
						SIMD dNX = load1(element.dN[gp][n][0]);
						SIMD dNY = load1(element.dN[gp][n][1]);
						SIMD dNZ = load1(element.dN[gp][n][2]);
						element.dND[gp][n][0] = inv0 * dNX + inv1 * dNY + inv2 * dNZ;
						element.dND[gp][n][1] = inv3 * dNX + inv4 * dNY + inv5 * dNZ;
						element.dND[gp][n][2] = inv6 * dNX + inv7 * dNY + inv8 * dNZ;
					}
				}
				//cylindrical
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD cooX =    element.gpcoords[gp][0];
					SIMD cooY =    element.gpcoords[gp][1];
					SIMD centerX = element.ecf.center[gp][0];
					SIMD centerY = element.ecf.center[gp][1];
					SIMD distanceX = cooX - centerX;
					SIMD distanceY = cooY - centerY;
					for (size_t s = 0; s < SIMD::size; ++s) {
						double rot = std::atan2(distanceY[s], distanceX[s]);
						element.cossin[gp][0][s] = 1;
						element.cossin[gp][1][s] = 1;
						element.cossin[gp][2][s] = std::cos(rot);
						element.cossin[gp][3][s] = 0;
						element.cossin[gp][4][s] = 0;
						element.cossin[gp][5][s] = std::sin(rot);
					}
				}
				//rotate cylindrical
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD cos0 = element.cossin[gp][0];
					SIMD cos1 = element.cossin[gp][1];
					SIMD cos2 = element.cossin[gp][2];
					SIMD sin0 = element.cossin[gp][3];
					SIMD sin1 = element.cossin[gp][4];
					SIMD sin2 = element.cossin[gp][5];

					SIMD t00 = cos1 * cos2;
					SIMD t01 = cos1 * sin2;
					SIMD t02 = -sin1;
					SIMD t10 = cos2 * sin0 * sin1 - cos0 * sin2;
					SIMD t11 = cos0 * cos2 + sin0 * sin1 * sin2;
					SIMD t12 = cos1 * sin0;
					SIMD t20 = sin0 * sin2 + cos0 * cos2 * sin1;
					SIMD t21 = cos0 * sin1 * sin2 - cos2 * sin0;
					SIMD t22 = cos0 * cos1;

					SIMD c00 = element.ecf.conductivity[gp][0];
					SIMD c01 = element.ecf.conductivity[gp][1], c11 = element.ecf.conductivity[gp][3];
					SIMD c02 = element.ecf.conductivity[gp][2], c12 = element.ecf.conductivity[gp][4], c22 = element.ecf.conductivity[gp][5];

					SIMD a = t00 * c00 + t10 * c01 + t20 * c02;
					SIMD b = t00 * c01 + t10 * c11 + t20 * c12;
					SIMD c = t00 * c02 + t10 * c12 + t20 * c22;
					element.conductivity[gp][0] = a * t00 + b * t10 + c * t20;
					element.conductivity[gp][1] = a * t01 + b * t11 + c * t21;
					element.conductivity[gp][2] = a * t02 + b * t12 + c * t22;

					a = t01 * c00 + t11 * c01 + t21 * c02;
					b = t01 * c01 + t11 * c11 + t21 * c12;
					c = t01 * c02 + t11 * c12 + t21 * c22;
					element.conductivity[gp][3] = a * t01 + b * t11 + c * t21;
					element.conductivity[gp][4] = a * t02 + b * t12 + c * t22;

					a = t02 * c00 + t12 * c01 + t22 * c02;
					b = t02 * c01 + t12 * c11 + t22 * c12;
					c = t02 * c02 + t12 * c12 + t22 * c22;
					element.conductivity[gp][5] = a * t02 + b * t12 + c * t22;
				}
				//stiffness symmetric non-isotropic
				double * __restrict__ out = stiffness.data;
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD scale = element.det[gp] * load1(element.w[gp]);
					SIMD c00 = element.conductivity[gp][0];
					SIMD c10 = element.conductivity[gp][1], c11 = element.conductivity[gp][3];
					SIMD c20 = element.conductivity[gp][2], c21 = element.conductivity[gp][4], c22 = element.conductivity[gp][5];
					for (size_t n = 0, i = 0; n < nodes; ++n) {
						SIMD nx = element.dND[gp][n][0];
						SIMD ny = element.dND[gp][n][1];
						SIMD nz = element.dND[gp][n][2];
						SIMD a = nx * c00 + ny * c10 + nz * c20;
						SIMD b = nx * c10 + ny * c11 + nz * c21;
						SIMD c = nx * c20 + ny * c21 + nz * c22;
						for (size_t m = n; m < nodes; ++m, ++i) {
							SIMD mx = element.dND[gp][m][0];
							SIMD my = element.dND[gp][m][1];
							SIMD mz = element.dND[gp][m][2];

							SIMD res = load(out + i * SIMD::size);
							res = res + scale * (a * mx + b * my + c * mz);
							store(out + i * SIMD::size, res);
						}
					}
				}
				stiffness+= SIMD::size;
			}
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xDEAD);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xFADE);
			}
			end = eslog::time();
		}

		if(cooToGP && rotateConductivity && isSpherical)
		{
			start = eslog::time();
			// std::cout<<"NON-CARTESIAN NON-ISOTROPIC YES ROTATION SYMMETRIC SPHERICAL -- SPHERICAL3D"<<std::endl;
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xFACE);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xCAFE);
			}
			esint chunks = elements / SIMD::size;
			for (esint c = 1; c < chunks; ++c) {
			// 	// cooToGP.simd(element);
				for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
					for (size_t n = 0; n < nodes; ++n) {
						for (size_t d = 0; d < ndim; ++d) {
							element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
						}
					}
				}
				for (size_t gp = 0; gp < gps; ++gp) {
					for (size_t d = 0; d < ndim; ++d) {
						element.gpcoords[gp][d] = zeros();
					}
					for (size_t n = 0; n < nodes; ++n) {
						for (size_t d = 0; d < ndim; ++d) {
							element.gpcoords[gp][d] = element.gpcoords[gp][d] + load1(element.N[gp][n]) * element.coords[n][d];
						}
					}
				}
				//integrate
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD jacobian0 = zeros(), jacobian1 = zeros(), jacobian2 = zeros(), jacobian3 = zeros(), jacobian4 = zeros();
					SIMD jacobian5 = zeros(), jacobian6 = zeros(), jacobian7 = zeros(), jacobian8 = zeros();

					for (size_t n = 0; n < nodes; ++n) {
						SIMD coordsX = element.coords[n][0];
						SIMD coordsY = element.coords[n][1];
						SIMD coordsZ = element.coords[n][2];
						SIMD dNX = load1(element.dN[gp][n][0]);
						SIMD dNY = load1(element.dN[gp][n][1]);
						SIMD dNZ = load1(element.dN[gp][n][2]);

						jacobian0 = jacobian0 + dNX * coordsX;
						jacobian1 = jacobian1 + dNX * coordsY;
						jacobian2 = jacobian2 + dNX * coordsZ;
						jacobian3 = jacobian3 + dNY * coordsX;
						jacobian4 = jacobian4 + dNY * coordsY;
						jacobian5 = jacobian5 + dNY * coordsZ;
						jacobian6 = jacobian6 + dNZ * coordsX;
						jacobian7 = jacobian7 + dNZ * coordsY;
						jacobian8 = jacobian8 + dNZ * coordsZ;
					}

					element.det[gp] =
							+ jacobian0 * jacobian4 * jacobian8
							+ jacobian1 * jacobian5 * jacobian6
							+ jacobian2 * jacobian3 * jacobian7
							- jacobian2 * jacobian4 * jacobian6
							- jacobian1 * jacobian3 * jacobian8
							- jacobian0 * jacobian5 * jacobian7;

					SIMD detJx = ones() / element.det[gp];
					SIMD inv0 = detJx * ( jacobian8 * jacobian4 - jacobian7 * jacobian5);
					SIMD inv1 = detJx * (-jacobian8 * jacobian1 + jacobian7 * jacobian2);
					SIMD inv2 = detJx * ( jacobian5 * jacobian1 - jacobian4 * jacobian2);
					SIMD inv3 = detJx * (-jacobian8 * jacobian3 + jacobian6 * jacobian5);
					SIMD inv4 = detJx * ( jacobian8 * jacobian0 - jacobian6 * jacobian2);
					SIMD inv5 = detJx * (-jacobian5 * jacobian0 + jacobian3 * jacobian2);
					SIMD inv6 = detJx * ( jacobian7 * jacobian3 - jacobian6 * jacobian4);
					SIMD inv7 = detJx * (-jacobian7 * jacobian0 + jacobian6 * jacobian1);
					SIMD inv8 = detJx * ( jacobian4 * jacobian0 - jacobian3 * jacobian1);

					for (size_t n = 0; n < nodes; ++n) {
						SIMD dNX = load1(element.dN[gp][n][0]);
						SIMD dNY = load1(element.dN[gp][n][1]);
						SIMD dNZ = load1(element.dN[gp][n][2]);
						element.dND[gp][n][0] = inv0 * dNX + inv1 * dNY + inv2 * dNZ;
						element.dND[gp][n][1] = inv3 * dNX + inv4 * dNY + inv5 * dNZ;
						element.dND[gp][n][2] = inv6 * dNX + inv7 * dNY + inv8 * dNZ;
					}
				}
				//spherical
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD x = element.gpcoords[gp][0] - element.ecf.center[gp][0];
					SIMD y = element.gpcoords[gp][1] - element.ecf.center[gp][1];
					SIMD z = element.gpcoords[gp][2] - element.ecf.center[gp][2];
					for (size_t s = 0; s < SIMD::size; ++s) {
						double azimut = std::atan2(y[s], x[s]);
						double r = std::sqrt(x[s] * x[s] + y[s] * y[s] + z[s] * z[s]);
						double elevation = r < 1e-15 ? 0 : std::atan2(std::sqrt(z[s] * z[s] + x[s] * x[s]), y[s]);
						element.cossin[gp][0][s] = 1;
						element.cossin[gp][1][s] = std::cos(elevation);
						element.cossin[gp][2][s] = std::cos(azimut);
						element.cossin[gp][3][s] = 0;
						element.cossin[gp][4][s] = std::sin(elevation);
						element.cossin[gp][5][s] = std::sin(azimut);
					}
				}
				//rotate spherical
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD cos0 = element.cossin[gp][0];
					SIMD cos1 = element.cossin[gp][1];
					SIMD cos2 = element.cossin[gp][2];
					SIMD sin0 = element.cossin[gp][3];
					SIMD sin1 = element.cossin[gp][4];
					SIMD sin2 = element.cossin[gp][5];

					SIMD t00 = cos1 * cos2;
					SIMD t01 = cos1 * sin2;
					SIMD t02 = -sin1;
					SIMD t10 = cos2 * sin0 * sin1 - cos0 * sin2;
					SIMD t11 = cos0 * cos2 + sin0 * sin1 * sin2;
					SIMD t12 = cos1 * sin0;
					SIMD t20 = sin0 * sin2 + cos0 * cos2 * sin1;
					SIMD t21 = cos0 * sin1 * sin2 - cos2 * sin0;
					SIMD t22 = cos0 * cos1;

					SIMD c00 = element.ecf.conductivity[gp][0];
					SIMD c01 = element.ecf.conductivity[gp][1], c11 = element.ecf.conductivity[gp][3];
					SIMD c02 = element.ecf.conductivity[gp][2], c12 = element.ecf.conductivity[gp][4], c22 = element.ecf.conductivity[gp][5];

					SIMD a = t00 * c00 + t10 * c01 + t20 * c02;
					SIMD b = t00 * c01 + t10 * c11 + t20 * c12;
					SIMD c = t00 * c02 + t10 * c12 + t20 * c22;
					element.conductivity[gp][0] = a * t00 + b * t10 + c * t20;
					element.conductivity[gp][1] = a * t01 + b * t11 + c * t21;
					element.conductivity[gp][2] = a * t02 + b * t12 + c * t22;

					a = t01 * c00 + t11 * c01 + t21 * c02;
					b = t01 * c01 + t11 * c11 + t21 * c12;
					c = t01 * c02 + t11 * c12 + t21 * c22;
					element.conductivity[gp][3] = a * t01 + b * t11 + c * t21;
					element.conductivity[gp][4] = a * t02 + b * t12 + c * t22;

					a = t02 * c00 + t12 * c01 + t22 * c02;
					b = t02 * c01 + t12 * c11 + t22 * c12;
					c = t02 * c02 + t12 * c12 + t22 * c22;
					element.conductivity[gp][5] = a * t02 + b * t12 + c * t22;
				}
				//stiffness symmetric non-isotropic
				double * __restrict__ out = stiffness.data;
				for (size_t gp = 0; gp < gps; ++gp) {
					SIMD scale = element.det[gp] * load1(element.w[gp]);
					SIMD c00 = element.conductivity[gp][0];
					SIMD c10 = element.conductivity[gp][1], c11 = element.conductivity[gp][3];
					SIMD c20 = element.conductivity[gp][2], c21 = element.conductivity[gp][4], c22 = element.conductivity[gp][5];
					for (size_t n = 0, i = 0; n < nodes; ++n) {
						SIMD nx = element.dND[gp][n][0];
						SIMD ny = element.dND[gp][n][1];
						SIMD nz = element.dND[gp][n][2];
						SIMD a = nx * c00 + ny * c10 + nz * c20;
						SIMD b = nx * c10 + ny * c11 + nz * c21;
						SIMD c = nx * c20 + ny * c21 + nz * c22;
						for (size_t m = n; m < nodes; ++m, ++i) {
							SIMD mx = element.dND[gp][m][0];
							SIMD my = element.dND[gp][m][1];
							SIMD mz = element.dND[gp][m][2];

							SIMD res = load(out + i * SIMD::size);
							res = res + scale * (a * mx + b * my + c * mz);
							store(out + i * SIMD::size, res);
						}
					}
				}
				stiffness+= SIMD::size;
			}
			if (action == ActionOperator::ASSEMBLE)
			{
				__SSC_MARK(0xDEAD);
			}
			if (action == ActionOperator::REASSEMBLE)
			{
				__SSC_MARK(0xFADE);
			}
			end = eslog::time();
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
	return measurements(initEnd - initStart, end - start);
}

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t ETYPE>
Assembler::measurements HeatTransfer::manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements,
	typename std::enable_if<
		ndim == 3 &&
		ETYPE == TransferElementType::ASYMMETRIC_ISOTROPIC, int>::type*
)
{
	return measurements();
}

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t ETYPE>
Assembler::measurements HeatTransfer::manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements,
	typename std::enable_if<
		ndim == 3 &&
		ETYPE == TransferElementType::ASYMMETRIC_GENERAL, int>::type*
)
{
	double initStart, initEnd;

	constexpr static double straightAngleRec = 1.0 / 180;

	initStart = eslog::time();
	eslog::info("       = LOOP TYPE                                                           MANUAL = \n");

	if (this->K == nullptr) {
		// std::cout<<"K IS NULL"<<std::endl;
		return loop<HeatTransferDataDescriptor, nodes, gps, ndim, edim, ETYPE>(action, ops, elements);
	}
	if ((action == ActionOperator::FILL) || (action == ActionOperator::SOLUTION))
	{
		return conditionsloop<DataDescriptor, nodes, gps, ndim, edim, ETYPE>(action, ops, interval, elements);
	}
	if (elements == 0) return measurements();

	typename DataDescriptor<nodes, gps, ndim, edim, ETYPE>::Element element;
	std::vector<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*> active; active.reserve(ops.size());

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if (elements > SIMD::size) {
				if ((*op)->isconst) {
					dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op)->simd(element);
				} else {
					active.push_back(dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op));
					active.back()->simd(element);
				}
			} else {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, ETYPE>*>(*op)->peel(element, elements);
			}
		}
	}

	auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;
	const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
	bool rotateConductivity = mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	bool isIsotropic = mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	// it is dirty hack just to be sure that compiler must assume both variants (currently settings.sigma = 0 and diffusion_split = false)
	bool constConductivity = !settings.diffusion_split;
	bool constRotation = settings.sigma == 0;
	if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
		if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			if (ndim == 2) {
				rotateConductivity &= mat->coordinate_system.rotation.z.isset;
			}
			if (ndim == 3) {
				rotateConductivity &= mat->coordinate_system.rotation.x.isset | mat->coordinate_system.rotation.y.isset | mat->coordinate_system.rotation.z.isset;
			}
		}
	}

	bool hasHeatSource = false;
	bool constHeatSource = true;
	auto heatSourceEval = configuration.heat_source.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (heatSourceEval != configuration.heat_source.end()) {
		hasHeatSource = true;
		constHeatSource = heatSourceEval->second.evaluator != nullptr;
	}

	bool hasAdvection = false;
	bool constAdvection = true;
	auto advectionEval = configuration.translation_motions.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (advectionEval != configuration.translation_motions.end()) {
		hasAdvection = true;
		constAdvection = advectionEval->second.x.evaluator != nullptr && advectionEval->second.y.evaluator != nullptr && advectionEval->second.z.evaluator != nullptr;
	}

	bool cooToGP = mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN;
	bool isSpherical = mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::SPHERICAL;
	bool computeK = action == ActionOperator::ASSEMBLE || action == ActionOperator::REASSEMBLE;
	bool computeGradient = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.gradient;
	bool computeFlux = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.flux;
	bool computeConductivity = computeK | computeFlux;
	bool getTemp = computeGradient || computeFlux;
	bool isfullMatrix = this->K->shape == Matrix_Shape::FULL;

	initEnd = eslog::time();

	OutputParameterIterator stiffness(this->elements.stiffness, interval);

	if((action == ActionOperator::VOID) || (action == ActionOperator::SOLUTION) || (action == ActionOperator::FILL))
	{
		std::cout<<"UNSUPPOERTED ACTION"<<std::endl;
		return measurements();
	}

	if( hasAdvection ||
		hasHeatSource ||
		getTemp)
	{
		std::cout<<"UNSUPPORTED OPERATOR"<<std::endl;
		return measurements();
	}

	if (!constRotation || !constConductivity)
	{
		std::cout<<"NON-CONST ROTATION OR CONDUCTIVITY"<<std::endl;
		return measurements();
	}

	// std::cout<<"ELEMENT TYPE IS "<<ETYPE<<std::endl;

	stiffness += SIMD::size;
	for (size_t s = 0; s < SIMD::size; ++s) {
		++procNodes;
	}

	// std::cout<<"WATCH OUT FOR THIS"<<std::endl;

	double start, end;
	if((action == ActionOperator::ASSEMBLE) || (action == ActionOperator::REASSEMBLE))
	{
		if(!cooToGP && rotateConductivity)
			{
				start = eslog::time();
				// std::cout<<"CARTESIAN NON-ISOTROPIC YES ROTATION NON-SYMETRIC -- ANISOTROPIC3D"<<std::endl;
				if (action == ActionOperator::ASSEMBLE)
				{
					__SSC_MARK(0xFACE);
				}
				if (action == ActionOperator::REASSEMBLE)
				{
					__SSC_MARK(0xCAFE);
				}
				esint chunks = elements / SIMD::size;
				for (esint c = 1; c < chunks; ++c) {
					//coo
					for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
						for (size_t n = 0; n < nodes; ++n) {
							for (size_t d = 0; d < ndim; ++d) {
								element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
							}
						}
					}
					//integrate
					for (size_t gp = 0; gp < gps; ++gp) {
						SIMD jacobian0 = zeros(), jacobian1 = zeros(), jacobian2 = zeros(), jacobian3 = zeros(), jacobian4 = zeros();
						SIMD jacobian5 = zeros(), jacobian6 = zeros(), jacobian7 = zeros(), jacobian8 = zeros();

						for (size_t n = 0; n < nodes; ++n) {
							SIMD coordsX = element.coords[n][0];
							SIMD coordsY = element.coords[n][1];
							SIMD coordsZ = element.coords[n][2];
							SIMD dNX = load1(element.dN[gp][n][0]);
							SIMD dNY = load1(element.dN[gp][n][1]);
							SIMD dNZ = load1(element.dN[gp][n][2]);

							jacobian0 = jacobian0 + dNX * coordsX;
							jacobian1 = jacobian1 + dNX * coordsY;
							jacobian2 = jacobian2 + dNX * coordsZ;
							jacobian3 = jacobian3 + dNY * coordsX;
							jacobian4 = jacobian4 + dNY * coordsY;
							jacobian5 = jacobian5 + dNY * coordsZ;
							jacobian6 = jacobian6 + dNZ * coordsX;
							jacobian7 = jacobian7 + dNZ * coordsY;
							jacobian8 = jacobian8 + dNZ * coordsZ;
						}

						element.det[gp] =
								+ jacobian0 * jacobian4 * jacobian8
								+ jacobian1 * jacobian5 * jacobian6
								+ jacobian2 * jacobian3 * jacobian7
								- jacobian2 * jacobian4 * jacobian6
								- jacobian1 * jacobian3 * jacobian8
								- jacobian0 * jacobian5 * jacobian7;

						SIMD detJx = ones() / element.det[gp];
						SIMD inv0 = detJx * ( jacobian8 * jacobian4 - jacobian7 * jacobian5);
						SIMD inv1 = detJx * (-jacobian8 * jacobian1 + jacobian7 * jacobian2);
						SIMD inv2 = detJx * ( jacobian5 * jacobian1 - jacobian4 * jacobian2);
						SIMD inv3 = detJx * (-jacobian8 * jacobian3 + jacobian6 * jacobian5);
						SIMD inv4 = detJx * ( jacobian8 * jacobian0 - jacobian6 * jacobian2);
						SIMD inv5 = detJx * (-jacobian5 * jacobian0 + jacobian3 * jacobian2);
						SIMD inv6 = detJx * ( jacobian7 * jacobian3 - jacobian6 * jacobian4);
						SIMD inv7 = detJx * (-jacobian7 * jacobian0 + jacobian6 * jacobian1);
						SIMD inv8 = detJx * ( jacobian4 * jacobian0 - jacobian3 * jacobian1);

						for (size_t n = 0; n < nodes; ++n) {
							SIMD dNX = load1(element.dN[gp][n][0]);
							SIMD dNY = load1(element.dN[gp][n][1]);
							SIMD dNZ = load1(element.dN[gp][n][2]);
							element.dND[gp][n][0] = inv0 * dNX + inv1 * dNY + inv2 * dNZ;
							element.dND[gp][n][1] = inv3 * dNX + inv4 * dNY + inv5 * dNZ;
							element.dND[gp][n][2] = inv6 * dNX + inv7 * dNY + inv8 * dNZ;
						}
					}
					// cartesian
					for (size_t gp = 0; gp < gps; ++gp) {
						SIMD angleX = element.ecf.center[gp][0];
						SIMD angleY = element.ecf.center[gp][1];
						SIMD angleZ = element.ecf.center[gp][2];
						for (size_t s = 0; s < SIMD::size; ++s) {
							element.cossin[gp][0][s] = std::cos(M_PI * angleX[s] * straightAngleRec);
							element.cossin[gp][1][s] = std::cos(M_PI * angleY[s] * straightAngleRec);
							element.cossin[gp][2][s] = std::cos(M_PI * angleZ[s] * straightAngleRec);
							element.cossin[gp][3][s] = std::sin(M_PI * angleX[s] * straightAngleRec);
							element.cossin[gp][4][s] = std::sin(M_PI * angleY[s] * straightAngleRec);
							element.cossin[gp][5][s] = std::sin(M_PI * angleZ[s] * straightAngleRec);
						}
					}
					//rotate cartesian
					for (size_t gp = 0; gp < gps; ++gp) {
						SIMD cos0 = element.cossin[gp][0];
						SIMD cos1 = element.cossin[gp][1];
						SIMD cos2 = element.cossin[gp][2];
						SIMD sin0 = element.cossin[gp][3];
						SIMD sin1 = element.cossin[gp][4];
						SIMD sin2 = element.cossin[gp][5];

						SIMD t00 = cos1 * cos2;
						SIMD t01 = cos1 * sin2;
						SIMD t02 = -sin1;
						SIMD t10 = cos2 * sin0 * sin1 - cos0 * sin2;
						SIMD t11 = cos0 * cos2 + sin0 * sin1 * sin2;
						SIMD t12 = cos1 * sin0;
						SIMD t20 = sin0 * sin2 + cos0 * cos2 * sin1;
						SIMD t21 = cos0 * sin1 * sin2 - cos2 * sin0;
						SIMD t22 = cos0 * cos1;

						SIMD origin0 = element.ecf.conductivity[gp][0];
						SIMD origin1 = element.ecf.conductivity[gp][1];
						SIMD origin2 = element.ecf.conductivity[gp][2];
						SIMD origin3 = element.ecf.conductivity[gp][3];
						SIMD origin4 = element.ecf.conductivity[gp][4];
						SIMD origin5 = element.ecf.conductivity[gp][5];
						SIMD origin6 = element.ecf.conductivity[gp][6];
						SIMD origin7 = element.ecf.conductivity[gp][7];
						SIMD origin8 = element.ecf.conductivity[gp][8];

						SIMD a = t00 * origin0 + t10 * origin3 + t20 * origin6;
						SIMD b = t00 * origin1 + t10 * origin4 + t20 * origin7;
						SIMD c = t00 * origin2 + t10 * origin5 + t20 * origin8;
						element.conductivity[gp][0] = a * t00 + b * t10 + c * t20;
						element.conductivity[gp][1] = a * t01 + b * t11 + c * t21;
						element.conductivity[gp][2] = a * t02 + b * t12 + c * t22;

						a = t01 * origin0 + t11 * origin3 + t21 * origin6;
						b = t01 * origin1 + t11 * origin4 + t21 * origin7;
						c = t01 * origin2 + t11 * origin5 + t21 * origin8;
						element.conductivity[gp][3] = a * t00 + b * t10 + c * t20;
						element.conductivity[gp][4] = a * t01 + b * t11 + c * t21;
						element.conductivity[gp][5] = a * t02 + b * t12 + c * t22;

						a = t02 * origin0 + t12 * origin3 + t22 * origin6;
						b = t02 * origin1 + t12 * origin4 + t22 * origin7;
						c = t02 * origin2 + t12 * origin5 + t22 * origin8;
						element.conductivity[gp][6] = a * t00 + b * t10 + c * t20;
						element.conductivity[gp][7] = a * t01 + b * t11 + c * t21;
						element.conductivity[gp][8] = a * t02 + b * t12 + c * t22;
					}
					//stiffness non-symmetric non-isotropic
					double * __restrict__ out = stiffness.data;
					for (size_t gp = 0; gp < gps; ++gp) {
						SIMD scale = element.det[gp] * load1(element.w[gp]);
						SIMD c00 = element.conductivity[gp][0], c01 = element.conductivity[gp][3], c02 = element.conductivity[gp][6];
						SIMD c10 = element.conductivity[gp][1], c11 = element.conductivity[gp][4], c12 = element.conductivity[gp][7];
						SIMD c20 = element.conductivity[gp][2], c21 = element.conductivity[gp][5], c22 = element.conductivity[gp][8];
						for (size_t n = 0; n < nodes; ++n) {
							SIMD nx = element.dND[gp][n][0];
							SIMD ny = element.dND[gp][n][1];
							SIMD nz = element.dND[gp][n][2];
							SIMD a = nx * c00 + ny * c01 + nz * c02;
							SIMD b = nx * c10 + ny * c11 + nz * c12;
							SIMD c = nx * c20 + ny * c21 + nz * c22;
							for (size_t m = 0; m < nodes; ++m) {
								SIMD mx = element.dND[gp][m][0];
								SIMD my = element.dND[gp][m][1];
								SIMD mz = element.dND[gp][m][2];
								SIMD res = load(out + (n * nodes + m) * SIMD::size);
								res = res + scale * (a * mx + b * my + c * mz);
								store(out + (n * nodes + m) * SIMD::size, res);
							}
						}
					}
					stiffness+= SIMD::size;

				}
				if (action == ActionOperator::ASSEMBLE)
				{
					__SSC_MARK(0xDEAD);
				}
				if (action == ActionOperator::REASSEMBLE)
				{
					__SSC_MARK(0xFADE);
				}
				end = eslog::time();
			}
	}

	// if(action == ActionOperator::REASSEMBLE)
	// {
	// 	start = eslog::time();
	// 	__SSC_MARK(0xCAFE);
	// 	esint chunks = elements / SIMD::size;
	// 	for (esint c = 1; c < chunks; ++c) {
	// 		if (cooToGP) {
	// 			cooAndGps.simd(element);
	// 		} else {
	// 			coo.simd(element);
	// 		}
	// 		integration.simd(element);
	// 		if (getTemp) {
	// 			temp.simd(element);
	// 		}
	// 		if (computeConductivity) {
	// 			if (!constConductivity) {
	// 				updateConductivity<DataDescriptor, nodes, gps, ndim, edim, ETYPE>()(element, mat);
	// 			}
	// 			if (rotateConductivity) {
	// 				if (!constRotation) {
	// 					updateRotation<DataDescriptor, nodes, gps, ndim, edim, ETYPE>()(element, mat);
	// 				}
	// 				rotation(element, mat);
	// 			}
	// 		}
	//
	// 		if (computeK) {
	// 			if (hasAdvection) {
	// 				if (!constAdvection) {
	// 					updateTM(element, advectionEval->second.x.evaluator);
	// 				}
	// 				advection.simd(element);
	// 			}
	// 			stiffness.simd(element);
	// 			if (hasHeatSource) {
	// 				if (!constHeatSource) {
	// 					updateHS(element, heatSourceEval->second.evaluator);
	// 				}
	// 				heatSource.simd(element);
	// 			}
	// 		}
	// 		if (action == ActionOperator::FILL) {
	// 			if (isfullMatrix) {
	// 				fullFiller.simd(element);
	// 			} else {
	// 				upperFiller.simd(element);
	// 			}
	// 			rhsFiller.simd(element);
	// 		}
	// 		if (computeGradient) {
	// 			gradient.simd(element);
	// 		}
	// 		if (computeFlux) {
	// 			flux.simd(element);
	// 		}
	// 	}
	// 	__SSC_MARK(0xFADE);
	// 	end = eslog::time();
	// }


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
	return measurements(initEnd - initStart, end - start);
}

template <int etype>
Assembler::measurements HeatTransfer::instantiateManual2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return manualloop<HeatTransferDataDescriptor, 3, HeatTransferGPC::TRIANGLE3, 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return manualloop<HeatTransferDataDescriptor, 6, HeatTransferGPC::TRIANGLE6, 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return manualloop<HeatTransferDataDescriptor, 4, HeatTransferGPC::SQUARE4  , 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return manualloop<HeatTransferDataDescriptor, 8, HeatTransferGPC::SQUARE8  , 2, 2, etype>(action, ops, interval, elements); break;
	default: return measurements();
	};
}

template <int etype>
Assembler::measurements HeatTransfer::instantiateManual3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TETRA4):    return manualloop<HeatTransferDataDescriptor,  4, HeatTransferGPC::TETRA4    , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::TETRA10):   return manualloop<HeatTransferDataDescriptor, 10, HeatTransferGPC::TETRA10   , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return manualloop<HeatTransferDataDescriptor,  5, HeatTransferGPC::PYRAMID5  , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): return manualloop<HeatTransferDataDescriptor, 13, HeatTransferGPC::PYRAMID13 , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA6):   return manualloop<HeatTransferDataDescriptor,  6, HeatTransferGPC::PRISMA6   , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA15):  return manualloop<HeatTransferDataDescriptor, 15, HeatTransferGPC::PRISMA15  , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA8):     return manualloop<HeatTransferDataDescriptor,  8, HeatTransferGPC::HEXA8     , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA20):    return manualloop<HeatTransferDataDescriptor, 20, HeatTransferGPC::HEXA20    , 3, 3, etype>(action, ops, interval, elements); break;
	default: return measurements();
	};
}

Assembler::measurements HeatTransfer::instantiateManual(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (etype) {
		// elements
		case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return instantiateManual2D<HeatTransferElementType::SYMMETRIC_ISOTROPIC >(action, code, ops, interval, elements);
		case HeatTransferElementType::SYMMETRIC_GENERAL   : return instantiateManual2D<HeatTransferElementType::SYMMETRIC_GENERAL   >(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return instantiateManual2D<HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_GENERAL  : return instantiateManual2D<HeatTransferElementType::ASYMMETRIC_GENERAL  >(action, code, ops, interval, elements);
		}
	case 3:
		switch (etype) {
		// elements
		case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return instantiateManual3D<HeatTransferElementType::SYMMETRIC_ISOTROPIC >(action, code, ops, interval, elements);
		case HeatTransferElementType::SYMMETRIC_GENERAL   : return instantiateManual3D<HeatTransferElementType::SYMMETRIC_GENERAL   >(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return instantiateManual3D<HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_GENERAL  : return instantiateManual3D<HeatTransferElementType::ASYMMETRIC_GENERAL  >(action, code, ops, interval, elements);
		}
	}
	return measurements();
}


}
