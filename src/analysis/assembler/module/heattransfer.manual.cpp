
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

					SIMD determinant = jacobian0 * jacobian3 - jacobian1 * jacobian2;

					SIMD detJx = ones() / determinant;
					SIMD inv0 =  detJx * jacobian3;
					SIMD inv1 = -detJx * jacobian1;
					SIMD inv2 = -detJx * jacobian2;
					SIMD inv3 =  detJx * jacobian0;
					
					SIMD scale = element.ecf.thickness[gp] * determinant * load1(element.w[gp]) * element.conductivity[gp];

					for (size_t n = 0; n < nodes; ++n) {
						SIMD dNX = load1(element.dN[gp][n][0]);
						SIMD dNY = load1(element.dN[gp][n][1]);
						element.dND[gp][n][0] = inv0 * dNX + inv1 * dNY;
						element.dND[gp][n][1] = inv2 * dNX + inv3 * dNY;
					}
					
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
	// eslog::info("       = LOOP TYPE                                                           MANUAL = \n");

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

					SIMD determinant = jacobian0 * jacobian3 - jacobian1 * jacobian2;

					SIMD detJx = ones() / determinant;
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

					SIMD c00 = element.conductivity[gp][0];
					SIMD c10 = element.conductivity[gp][1], c11 = element.conductivity[gp][2];
					SIMD scale = element.ecf.thickness[gp] * determinant * load1(element.w[gp]);
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
			SIMD angle = element.ecf.center[0][0];
			SIMD cos;
			SIMD sin;
			for (size_t s = 0; s < SIMD::size; ++s) {
				cos[s] = std::cos(M_PI * angle[s] * straightAngleRec);
				sin[s] = std::sin(M_PI * angle[s] * straightAngleRec);
			}//TODO: only 1 gp of 1 el
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

					SIMD determinant = jacobian0 * jacobian3 - jacobian1 * jacobian2;

					SIMD detJx = ones() / determinant;
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


					SIMD ic00 = element.ecf.conductivity[gp][0];
					SIMD ic10 = element.ecf.conductivity[gp][1], ic11 = element.ecf.conductivity[gp][2];

					SIMD c00;
					SIMD c10, c11;

					c00 = (cos * ic00 - sin * ic10) * cos - (cos * ic10 - sin * ic11) * sin;
					c10 = (cos * ic00 - sin * ic10) * sin + (cos * ic10 - sin * ic11) * cos;
					c11 = (sin * ic00 + cos * ic10) * sin + (sin * ic10 + cos * ic11) * cos;


					SIMD scale = element.ecf.thickness[gp] * determinant * load1(element.w[gp]);
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
				esint chunks = elements / SIMD::size;
				std::vector<double> &storage = cossin_conditions[interval];
				double* iterator;
				storage.resize(2 * elements * gps * SIMD::size);
				iterator = storage.data();

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

						SIMD determinant = jacobian0 * jacobian3 - jacobian1 * jacobian2;

						SIMD detJx = ones() / determinant;
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

						SIMD cooX =    element.gpcoords[gp][0];
						SIMD cooY =    element.gpcoords[gp][1];
						SIMD centerX = element.ecf.center[gp][0];
						SIMD centerY = element.ecf.center[gp][1];
						SIMD distanceX = cooX - centerX;
						SIMD distanceY = cooY - centerY;
						SIMD cossin[2];
						for (size_t s = 0; s < SIMD::size; ++s) {
							double rot = std::atan2(distanceY[s], distanceX[s]);
							cossin[0][s] = std::cos(rot);
							cossin[1][s] = std::sin(rot);
						}

						memcpy(iterator, cossin, 2 * SIMD::size * sizeof(double));
						iterator += 2 * SIMD::size;

						SIMD ic00 = element.ecf.conductivity[gp][0];
						SIMD ic10 = element.ecf.conductivity[gp][1], ic11 = element.ecf.conductivity[gp][2];

						SIMD c00;
						SIMD c10, c11;

						c00 = (cossin[0] * ic00 - cossin[1] * ic10) * cossin[0] - (cossin[0] * ic10 - cossin[1] * ic11) * cossin[1];
						c10 = (cossin[0] * ic00 - cossin[1] * ic10) * cossin[1] + (cossin[0] * ic10 - cossin[1] * ic11) * cossin[0];
						c11 = (cossin[1] * ic00 + cossin[0] * ic10) * cossin[1] + (cossin[1] * ic10 + cossin[0] * ic11) * cossin[0];

						SIMD scale = element.ecf.thickness[gp] * determinant * load1(element.w[gp]);
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
	// eslog::info("       = LOOP TYPE                                                           MANUAL = \n");

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
			SIMD angle = element.ecf.center[0][0];
			SIMD cos;
			SIMD sin;
			for (size_t s = 0; s < SIMD::size; ++s) {
				cos[s] = std::cos(M_PI * angle[s] * straightAngleRec);
				sin[s] = std::sin(M_PI * angle[s] * straightAngleRec);
			}
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

					SIMD determinant = jacobian0 * jacobian3 - jacobian1 * jacobian2;

					SIMD detJx = ones() / determinant;
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

					SIMD origin0 = element.ecf.conductivity[gp][0];
					SIMD origin1 = element.ecf.conductivity[gp][1];
					SIMD origin2 = element.ecf.conductivity[gp][2];
					SIMD origin3 = element.ecf.conductivity[gp][3];
					

					SIMD c00 = (cos * origin0 - sin * origin2) * cos - (cos * origin1 - sin * origin3) * sin;
					SIMD c10 = (cos * origin0 - sin * origin2) * sin + (cos * origin1 - sin * origin3) * cos;
					SIMD c01 = (sin * origin0 + cos * origin2) * cos - (sin * origin1 + cos * origin3) * sin;
					SIMD c11 = (sin * origin0 + cos * origin2) * sin + (sin * origin1 + cos * origin3) * cos;


					SIMD scale = element.ecf.thickness[gp] * determinant * load1(element.w[gp]);
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
	// eslog::info("       = LOOP TYPE                                                           MANUAL = \n");

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
				
				double * __restrict__ out = stiffness.data;
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

					SIMD determinant =
							+ jacobian0 * jacobian4 * jacobian8
							+ jacobian1 * jacobian5 * jacobian6
							+ jacobian2 * jacobian3 * jacobian7
							- jacobian2 * jacobian4 * jacobian6
							- jacobian1 * jacobian3 * jacobian8
							- jacobian0 * jacobian5 * jacobian7;

					SIMD detJx = ones() / determinant;
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

					SIMD scale = determinant * load1(element.w[gp]) * element.conductivity[gp];
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
	// eslog::info("       = LOOP TYPE                                                           MANUAL = \n");

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
				double * __restrict__ out = stiffness.data;
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

					SIMD determinant =
							+ jacobian0 * jacobian4 * jacobian8
							+ jacobian1 * jacobian5 * jacobian6
							+ jacobian2 * jacobian3 * jacobian7
							- jacobian2 * jacobian4 * jacobian6
							- jacobian1 * jacobian3 * jacobian8
							- jacobian0 * jacobian5 * jacobian7;

					SIMD detJx = ones() / determinant;
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

					//stiffness symmetric non-isotropic
					SIMD scale = determinant * load1(element.w[gp]);
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
			SIMD angleX = element.ecf.center[0][0];
			SIMD angleY = element.ecf.center[0][1];
			SIMD angleZ = element.ecf.center[0][2];
			SIMD cos0;
			SIMD cos1;
			SIMD cos2;
			SIMD sin0;
			SIMD sin1;
			SIMD sin2;

			for (size_t s = 0; s < SIMD::size; ++s) {
				cos0[s] = std::cos(M_PI * angleX[s] * straightAngleRec);
				cos1[s] = std::cos(M_PI * angleY[s] * straightAngleRec);
				cos2[s] = std::cos(M_PI * angleZ[s] * straightAngleRec);
				sin0[s] = std::sin(M_PI * angleX[s] * straightAngleRec);
				sin1[s] = std::sin(M_PI * angleY[s] * straightAngleRec);
				sin2[s] = std::sin(M_PI * angleZ[s] * straightAngleRec);
			}//TODO only 1ce for 1 gp and 1 el
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
				double * __restrict__ out = stiffness.data;
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

					SIMD determinant =
							+ jacobian0 * jacobian4 * jacobian8
							+ jacobian1 * jacobian5 * jacobian6
							+ jacobian2 * jacobian3 * jacobian7
							- jacobian2 * jacobian4 * jacobian6
							- jacobian1 * jacobian3 * jacobian8
							- jacobian0 * jacobian5 * jacobian7;

					SIMD detJx = ones() / determinant;
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

					SIMD t00 = cos1 * cos2;
					SIMD t01 = cos1 * sin2;
					SIMD t02 = -sin1;
					SIMD t10 = cos2 * sin0 * sin1 - cos0 * sin2;
					SIMD t11 = cos0 * cos2 + sin0 * sin1 * sin2;
					SIMD t12 = cos1 * sin0;
					SIMD t20 = sin0 * sin2 + cos0 * cos2 * sin1;
					SIMD t21 = cos0 * sin1 * sin2 - cos2 * sin0;
					SIMD t22 = cos0 * cos1;

					SIMD ic00 = element.ecf.conductivity[gp][0];
					SIMD ic10 = element.ecf.conductivity[gp][1], ic11 = element.ecf.conductivity[gp][3];
					SIMD ic20 = element.ecf.conductivity[gp][2], ic21 = element.ecf.conductivity[gp][4], ic22 = element.ecf.conductivity[gp][5];

					SIMD c00;
					SIMD c10, c11;
					SIMD c20, c21, c22;

					SIMD a = t00 * ic00 + t10 * ic10 + t20 * ic20;
					SIMD b = t00 * ic10 + t10 * ic11 + t20 * ic21;
					SIMD c = t00 * ic20 + t10 * ic21 + t20 * ic22;
					c00 = a * t00 + b * t10 + c * t20;
					c10 = a * t01 + b * t11 + c * t21;
					c20 = a * t02 + b * t12 + c * t22;

					a = t01 * ic00 + t11 * ic10 + t21 * ic20;
					b = t01 * ic10 + t11 * ic11 + t21 * ic21;
					c = t01 * ic20 + t11 * ic21 + t21 * ic22;
					c11 = a * t01 + b * t11 + c * t21;
					c21 = a * t02 + b * t12 + c * t22;

					a = t02 * ic00 + t12 * ic10 + t22 * ic20;
					b = t02 * ic10 + t12 * ic11 + t22 * ic21;
					c = t02 * ic20 + t12 * ic21 + t22 * ic22;
					c22 = a * t02 + b * t12 + c * t22;

					SIMD scale = determinant * load1(element.w[gp]);
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
			esint chunks = elements / SIMD::size;
			std::vector<double> &storage = cossin_conditions[interval];
			double* iterator;
			storage.resize(6 * elements * gps * SIMD::size);
			iterator = storage.data();
			
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
				double * __restrict__ out = stiffness.data;
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

					SIMD determinant =
							+ jacobian0 * jacobian4 * jacobian8
							+ jacobian1 * jacobian5 * jacobian6
							+ jacobian2 * jacobian3 * jacobian7
							- jacobian2 * jacobian4 * jacobian6
							- jacobian1 * jacobian3 * jacobian8
							- jacobian0 * jacobian5 * jacobian7;

					SIMD detJx = ones() / determinant;
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
					
					//cylindrical
					SIMD cooX =    element.gpcoords[gp][0];
					SIMD cooY =    element.gpcoords[gp][1];
					SIMD centerX = element.ecf.center[gp][0];
					SIMD centerY = element.ecf.center[gp][1];
					SIMD distanceX = cooX - centerX;
					SIMD distanceY = cooY - centerY;
					SIMD cossin[6];
					for (size_t s = 0; s < SIMD::size; ++s) {
						double rot = std::atan2(distanceY[s], distanceX[s]);
						cossin[0][s] = 1;
						cossin[1][s] = 1;
						cossin[2][s] = std::cos(rot);
						cossin[3][s] = 0;
						cossin[4][s] = 0;
						cossin[5][s] = std::sin(rot);
					}
					
					memcpy(iterator, cossin, 6 * SIMD::size * sizeof(double));
					iterator += 6 * SIMD::size;

					//rotate cylindrical
					SIMD t00 = cossin[1] * cossin[2];
					SIMD t01 = cossin[1] * cossin[5];
					SIMD t02 = -cossin[4];
					SIMD t10 = cossin[2] * cossin[3] * cossin[4] - cossin[0] * cossin[5];
					SIMD t11 = cossin[0] * cossin[2] + cossin[3] * cossin[4] * cossin[5];
					SIMD t12 = cossin[1] * cossin[3];
					SIMD t20 = cossin[3] * cossin[5] + cossin[0] * cossin[2] * cossin[4];
					SIMD t21 = cossin[0] * cossin[4] * cossin[5] - cossin[2] * cossin[3];
					SIMD t22 = cossin[0] * cossin[1];

					SIMD ic00 = element.ecf.conductivity[gp][0];
					SIMD ic10 = element.ecf.conductivity[gp][1], ic11 = element.ecf.conductivity[gp][3];
					SIMD ic20 = element.ecf.conductivity[gp][2], ic21 = element.ecf.conductivity[gp][4], ic22 = element.ecf.conductivity[gp][5];

					SIMD c00;
					SIMD c10, c11;
					SIMD c20, c21, c22;

					SIMD a = t00 * ic00 + t10 * ic10 + t20 * ic20;
					SIMD b = t00 * ic10 + t10 * ic11 + t20 * ic21;
					SIMD c = t00 * ic20 + t10 * ic21 + t20 * ic22;
					c00 = a * t00 + b * t10 + c * t20;
					c10 = a * t01 + b * t11 + c * t21;
					c20 = a * t02 + b * t12 + c * t22;

					a = t01 * ic00 + t11 * ic10 + t21 * ic20;
					b = t01 * ic10 + t11 * ic11 + t21 * ic21;
					c = t01 * ic20 + t11 * ic21 + t21 * ic22;
					c11 = a * t01 + b * t11 + c * t21;
					c21 = a * t02 + b * t12 + c * t22;

					a = t02 * ic00 + t12 * ic10 + t22 * ic20;
					b = t02 * ic10 + t12 * ic11 + t22 * ic21;
					c = t02 * ic20 + t12 * ic21 + t22 * ic22;
					c22 = a * t02 + b * t12 + c * t22;

					//stiffness symmetric non-isotropic
					SIMD scale = determinant * load1(element.w[gp]);

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
			esint chunks = elements / SIMD::size;
			std::vector<double> &storage = cossin_conditions[interval];
			double* iterator;
			storage.resize(6 * elements * gps * SIMD::size);
			iterator = storage.data();

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
				double * __restrict__ out = stiffness.data;
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

					SIMD determinant =
							+ jacobian0 * jacobian4 * jacobian8
							+ jacobian1 * jacobian5 * jacobian6
							+ jacobian2 * jacobian3 * jacobian7
							- jacobian2 * jacobian4 * jacobian6
							- jacobian1 * jacobian3 * jacobian8
							- jacobian0 * jacobian5 * jacobian7;

					SIMD detJx = ones() / determinant;
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
					//spherical
					SIMD x = element.gpcoords[gp][0] - element.ecf.center[gp][0];
					SIMD y = element.gpcoords[gp][1] - element.ecf.center[gp][1];
					SIMD z = element.gpcoords[gp][2] - element.ecf.center[gp][2];
					SIMD cossin[6];

					for (size_t s = 0; s < SIMD::size; ++s) {
						double azimut = std::atan2(y[s], x[s]);
						double r = std::sqrt(x[s] * x[s] + y[s] * y[s] + z[s] * z[s]);
						double elevation = r < 1e-15 ? 0 : std::atan2(std::sqrt(z[s] * z[s] + x[s] * x[s]), y[s]);
						cossin[0][s] = 1;
						cossin[1][s] = std::cos(elevation);
						cossin[2][s] = std::cos(azimut);
						cossin[3][s] = 0;
						cossin[4][s] = std::sin(elevation);
						cossin[5][s] = std::sin(azimut);
					}
					
					memcpy(iterator, cossin, 6 * SIMD::size * sizeof(double));
					iterator += 6 * SIMD::size;
					
					//rotate spherical

					SIMD t00 = cossin[1] * cossin[2];
					SIMD t01 = cossin[1] * cossin[5];
					SIMD t02 = -cossin[4];
					SIMD t10 = cossin[2] * cossin[3] * cossin[4] - cossin[0] * cossin[5];
					SIMD t11 = cossin[0] * cossin[2] + cossin[3] * cossin[4] * cossin[5];
					SIMD t12 = cossin[1] * cossin[3];
					SIMD t20 = cossin[3] * cossin[5] + cossin[0] * cossin[2] * cossin[4];
					SIMD t21 = cossin[0] * cossin[4] * cossin[5] - cossin[2] * cossin[3];
					SIMD t22 = cossin[0] * cossin[1];

					SIMD ic00 = element.ecf.conductivity[gp][0];
					SIMD ic10 = element.ecf.conductivity[gp][1], ic11 = element.ecf.conductivity[gp][3];
					SIMD ic20 = element.ecf.conductivity[gp][2], ic21 = element.ecf.conductivity[gp][4], ic22 = element.ecf.conductivity[gp][5];

					SIMD c00;
					SIMD c10, c11;
					SIMD c20, c21, c22;

					SIMD a = t00 * ic00 + t10 * ic10 + t20 * ic20;
					SIMD b = t00 * ic10 + t10 * ic11 + t20 * ic21;
					SIMD c = t00 * ic20 + t10 * ic21 + t20 * ic22;
					c00 = a * t00 + b * t10 + c * t20;
					c10 = a * t01 + b * t11 + c * t21;
					c20 = a * t02 + b * t12 + c * t22;

					a = t01 * ic00 + t11 * ic10 + t21 * ic20;
					b = t01 * ic10 + t11 * ic11 + t21 * ic21;
					c = t01 * ic20 + t11 * ic21 + t21 * ic22;
					c11 = a * t01 + b * t11 + c * t21;
					c21 = a * t02 + b * t12 + c * t22;

					a = t02 * ic00 + t12 * ic10 + t22 * ic20;
					b = t02 * ic10 + t12 * ic11 + t22 * ic21;
					c = t02 * ic20 + t12 * ic21 + t22 * ic22;
					c22 = a * t02 + b * t12 + c * t22;
					
					//stiffness symmetric non-isotropic
					SIMD scale = determinant * load1(element.w[gp]);

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
	// eslog::info("       = LOOP TYPE                                                           MANUAL = \n");

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
			// cartesian
			SIMD angleX = element.ecf.center[0][0];
			SIMD angleY = element.ecf.center[0][1];
			SIMD angleZ = element.ecf.center[0][2];
			SIMD cos0;
			SIMD cos1;
			SIMD cos2;
			SIMD sin0;
			SIMD sin1;
			SIMD sin2;
			for (size_t s = 0; s < SIMD::size; ++s) {
				cos0[s] = std::cos(M_PI * angleX[s] * straightAngleRec);
				cos1[s] = std::cos(M_PI * angleY[s] * straightAngleRec);
				cos2[s] = std::cos(M_PI * angleZ[s] * straightAngleRec);
				sin0[s] = std::sin(M_PI * angleX[s] * straightAngleRec);
				sin1[s] = std::sin(M_PI * angleY[s] * straightAngleRec);
				sin2[s] = std::sin(M_PI * angleZ[s] * straightAngleRec);
			}
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
				double * __restrict__ out = stiffness.data;
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

					SIMD determinant =
							+ jacobian0 * jacobian4 * jacobian8
							+ jacobian1 * jacobian5 * jacobian6
							+ jacobian2 * jacobian3 * jacobian7
							- jacobian2 * jacobian4 * jacobian6
							- jacobian1 * jacobian3 * jacobian8
							- jacobian0 * jacobian5 * jacobian7;

					SIMD detJx = ones() / determinant;
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
					SIMD c00 = a * t00 + b * t10 + c * t20;
					SIMD c10 = a * t01 + b * t11 + c * t21;
					SIMD c20 = a * t02 + b * t12 + c * t22;

					a = t01 * origin0 + t11 * origin3 + t21 * origin6;
					b = t01 * origin1 + t11 * origin4 + t21 * origin7;
					c = t01 * origin2 + t11 * origin5 + t21 * origin8;
					SIMD c01 = a * t00 + b * t10 + c * t20;
					SIMD c11 = a * t01 + b * t11 + c * t21;
					SIMD c21 = a * t02 + b * t12 + c * t22;

					a = t02 * origin0 + t12 * origin3 + t22 * origin6;
					b = t02 * origin1 + t12 * origin4 + t22 * origin7;
					c = t02 * origin2 + t12 * origin5 + t22 * origin8;
					SIMD c02 = a * t00 + b * t10 + c * t20;
					SIMD c12 = a * t01 + b * t11 + c * t21;
					SIMD c22 = a * t02 + b * t12 + c * t22;

					//stiffness non-symmetric non-isotropic
					SIMD scale = determinant * load1(element.w[gp]);

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
