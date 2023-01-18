
#include "expression.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/assembler/module/acoustic.h"
#include "analysis/assembler/module/structuralmechanics.h"
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

template <typename Module, template<size_t, size_t> class Operator>
void _fromExpression(Module &module, ParameterData &parameter, ExternalElementValue &value)
{
	if (std::all_of(value.evaluator.begin(), value.evaluator.end(), [] (const Evaluator *ev) { return ev == NULL; })) {
		return;
	}

	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		for (int d = 0; d < value.dimension; ++d) {
			if (value.evaluator[i * value.dimension + d]) {
				parameter.update[i] = 1;
				for (size_t p = 0; p < value.evaluator[i * value.dimension + d]->params.general.size(); ++p) {
					module.controller.addInput(i, parameter, value.evaluator[i * value.dimension + d]->params.general[p].variable);
				}
			}
		}
	}

	module.controller.prepare(parameter);

	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		for (int d = 0; d < value.dimension; ++d) {
			if (value.evaluator[i * value.dimension + d]) {
				module.elementOps[i].emplace_back(instantiate<typename Module::NGP, Operator>(i, module.controller, parameter, value.evaluator[i * value.dimension + d], d, value.dimension));
			}
		}
	}
}


void fromExpression(HeatTransfer &module, ParameterData &parameter, ExternalElementNodesValue &value)
{
	_fromExpression<HeatTransfer, ExpressionsToNodes>(module, parameter, value);
}

void fromExpression(HeatTransfer &module, ParameterData &parameter, ExternalElementGPsValue &value)
{
	_fromExpression<HeatTransfer, ExpressionsToGPs>(module, parameter, value);
}

void fromExpression(Acoustic &module, ParameterData &parameter, ExternalElementGPsValue &value)
{
	_fromExpression<Acoustic, ExpressionsToGPs>(module, parameter, value);
}

void fromExpression(StructuralMechanics &module, ParameterData &parameter, ExternalElementNodesValue &value)
{
	_fromExpression<StructuralMechanics, ExpressionsToNodes>(module, parameter, value);
}

void fromExpression(StructuralMechanics &module, ParameterData &parameter, ExternalElementGPsValue &value)
{
	_fromExpression<StructuralMechanics, ExpressionsToGPs>(module, parameter, value);
}


template <typename Module>
void _fromExpression(Module &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values)
{
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		for (int d = 0; d < values.dimension; ++d) {
			if (values.evaluator[r * values.dimension + d]) {
				for (size_t p = 0; p < values.evaluator[r * values.dimension + d]->params.general.size(); ++p) {
					module.controller.addInput(parameter.regions[r], values.evaluator[r * values.dimension + d]->params.general[p].variable);
				}
			}
		}

		if (std::any_of(values.evaluator.begin() + r * values.dimension, values.evaluator.begin() + r * values.dimension + values.dimension, [] (const Evaluator *ev) { return ev != NULL; })) {
			module.controller.prepare(parameter.regions[r]);
			for (int d = 0; d < values.dimension; ++d) {
				if (values.evaluator[r * values.dimension + d]) {
					std::fill(parameter.regions[r].update.begin(), parameter.regions[r].update.end(), 1);
					if (info::mesh->boundaryRegions[r]->dimension) {
						for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
							module.boundaryOps[r][i].emplace_back(instantiate<typename Module::NGP, ExpressionsToGPs>(r, i, module.controller, parameter.regions[r], values.evaluator[r * values.dimension + d], d, values.dimension));
						}
					} else {
						for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
							module.boundaryOps[r][t].emplace_back(instantiate<typename Module::NGP, ExpressionsToNodes>(r, t, module.controller, parameter.regions[r], values.evaluator[r * values.dimension + d], d, values.dimension));
						}
					}
				}
			}
		}
	}
}

void fromExpression(HeatTransfer &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values)
{
	_fromExpression(module, parameter, values);
}

void fromExpression(Acoustic &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values)
{
	_fromExpression(module, parameter, values);
}

void fromExpression(StructuralMechanics &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values)
{
	_fromExpression(module, parameter, values);
}

}
