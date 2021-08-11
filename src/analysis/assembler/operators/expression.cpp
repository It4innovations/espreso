
#include "expression.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/assembler/module/acoustic.h"
#include "basis/evaluator/evaluator.h"
#include "basis/expression/variable.h"
#include "basis/utilities/parser.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

#include <algorithm>
#include <memory>

namespace espreso {

template <template<size_t, size_t> class Operator>
void _fromExpression(AX_HeatTransfer &module, ParameterData &parameter, ExternalElementValue &value)
{
	if (std::all_of(value.evaluator.begin(), value.evaluator.end(), [] (const Evaluator *ev) { return ev == NULL; })) {
		return;
	}

	parameter.resize();
	module.addParameter(parameter);

	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		for (int d = 0; d < value.dimension && value.evaluator[i * value.dimension + d]; ++d) {
			module.elementOps[i].emplace_back(instantiate<AX_HeatTransfer::NGP, Operator>(i, parameter, value.evaluator[i * value.dimension + d], d, value.dimension));
		}
	}
}

void fromExpression(AX_HeatTransfer &module, ParameterData &parameter, ExternalElementNodesValue &value)
{
	_fromExpression<ExpressionsToNodes>(module, parameter, value);
}

void fromExpression(AX_HeatTransfer &module, ParameterData &parameter, ExternalElementGPsValue &value)
{
	_fromExpression<ExpressionsToGPs>(module, parameter, value);
}

template <template<size_t, size_t> typename Operator>
void _fromExpression(AX_HeatTransfer &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values)
{
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		for (int d = 0; d < values.dimension && values.evaluator[r * values.dimension + d]; ++d) {
			// add input parameters
		}

		if (std::any_of(values.evaluator.begin() + r * values.dimension, values.evaluator.begin() + r * values.dimension + values.dimension, [] (const Evaluator *ev) { return ev != NULL; })) {
			parameter.regions[r].resize();
			module.addParameter(parameter.regions[r]);
			for (int d = 0; d < values.dimension && values.evaluator[r * values.dimension + d]; ++d) {
				for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
					module.boundaryOps[r][i].emplace_back(instantiate<AX_HeatTransfer::NGP, Operator>(r, i, parameter.regions[r], values.evaluator[r * values.dimension + d], d, values.dimension));
				}
			}
		}
	}
}

void fromExpression(AX_HeatTransfer &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values)
{
	_fromExpression<ExpressionsToGPs>(module, parameter, values);
}

template <template<size_t, size_t> typename Operator>
void _evaluateFromExpression(AX_HeatTransfer &module, ParameterData &parameter, ExternalElementValue &value)
{
	if (std::all_of(value.evaluator.begin(), value.evaluator.end(), [] (const Evaluator *ev) { return ev == NULL; })) {
		return;
	}

	parameter.resize();
	module.addParameter(parameter);

	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		for (int d = 0; d < value.dimension && value.evaluator[i * value.dimension + d]; ++d) {
//			if (value.evaluator[i * value.dimension + d]->)
			std::unique_ptr<ActionOperator> op(instantiate<AX_HeatTransfer::NGP, Operator>(i, parameter, value.evaluator[i * value.dimension + d], d, value.dimension));
			size_t elementsInInterval = info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;

			for (size_t element = 0; element < elementsInInterval; ++element) {
				if(element == 0 || !op->isconst) {
					(*op)();
					++(*op);
				}
			}
		}
	}
}

void evaluateFromExpression(AX_HeatTransfer &module, ParameterData &parameter, ExternalElementNodesValue &value)
{
	_evaluateFromExpression<ExpressionsToNodes>(module, parameter, value);
}

void evaluateFromExpression(AX_HeatTransfer &module, ParameterData &parameter, ExternalElementGPsValue &value)
{
	_evaluateFromExpression<ExpressionsToGPs>(module, parameter, value);
}

void fromExpression(AX_Acoustic &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values)
{
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		for (int d = 0; d < values.dimension && values.evaluator[r * values.dimension + d]; ++d) {
			// add input parameters
		}

		if (std::any_of(values.evaluator.begin() + r * values.dimension, values.evaluator.begin() + r * values.dimension + values.dimension, [] (const Evaluator *ev) { return ev != NULL; })) {
			parameter.regions[r].resize();
			module.addParameter(parameter.regions[r]);
			for (int d = 0; d < values.dimension && values.evaluator[r * values.dimension + d]; ++d) {
				for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
					module.boundaryOps[r][i].emplace_back(instantiate<AX_HeatTransfer::NGP, ExpressionsToGPs>(r, i, parameter.regions[r], values.evaluator[r * values.dimension + d], d, values.dimension));
				}
			}
		}
	}
}

}
