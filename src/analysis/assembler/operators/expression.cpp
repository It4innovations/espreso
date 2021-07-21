
#include "expression.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/module/heattransfer.h"
#include "basis/evaluator/evaluator.h"
#include "basis/utilities/parser.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

#include <algorithm>

namespace espreso {

void fromExpression(AX_HeatTransfer &module, ParameterData &parameter, ExternalElementValue &value)
{
	if (std::all_of(value.evaluator.begin(), value.evaluator.end(), [] (const Evaluator *ev) { return ev == NULL; })) {
		return;
	}

	for (size_t i = 0; i < parameter.isconst.size(); ++i) {
		for (int d = 0; d < value.dimension && value.evaluator[i * value.dimension + d]; ++d) {
			for (size_t p = 0; p < value.evaluator[i * value.dimension + d]->variables.size(); ++p) {
				parameter.isconst[i] = false; // in the case of TIME it is possible to keep value constant
				if (StringCompare::caseInsensitiveEq("INITIAL_TEMPERATURE", value.evaluator[i * value.dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						value.evaluator[i * value.dimension + d]->params.general.push_back({ (module.temp.initial.gp.data->begin() + i)->data(), 0, module.temp.initial.gp.isconst[i] ? 0 : 1 });
						parameter.addInput(module.temp.initial.gp);
					} else {
						value.evaluator[i * value.dimension + d]->params.general.push_back({ (module.temp.initial.node.data->begin() + i)->data(), 0, module.temp.initial.node.isconst[i] ? 0 : 1 });
						parameter.addInput(module.temp.initial.node);
					}
				}
				if (StringCompare::caseInsensitiveEq("TEMPERATURE", value.evaluator[i * value.dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						value.evaluator[i * value.dimension + d]->params.general.push_back({ (module.temp.gp.data->begin() + i)->data(), 0, module.temp.gp.isconst[i] ? 0 : 1 });
						parameter.addInput(module.temp.gp);
					} else {
						value.evaluator[i * value.dimension + d]->params.general.push_back({ (module.temp.node.data->begin() + i)->data(), 0, module.temp.node.isconst[i] ? 0 : 1 });
						parameter.addInput(module.temp.node);
					}
				}
				bool insertGpCoords = false, insertNodeCoords = false;
				if (StringCompare::caseInsensitiveEq("COORDINATE_X", value.evaluator[i * value.dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						value.evaluator[i * value.dimension + d]->params.general.push_back({ (module.coords.gp.data->begin() + i)->data(), 0, module.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
						insertGpCoords = true;
					} else {
						value.evaluator[i * value.dimension + d]->params.general.push_back({ (module.coords.node.data->begin() + i)->data(), 0, module.coords.node.isconst[i] ? 0 : info::mesh->dimension });
						insertNodeCoords = true;
					}
				}
				if (StringCompare::caseInsensitiveEq("COORDINATE_Y", value.evaluator[i * value.dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						value.evaluator[i * value.dimension + d]->params.general.push_back({ (module.coords.gp.data->begin() + i)->data(), 1, module.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
						insertGpCoords = true;
					} else {
						value.evaluator[i * value.dimension + d]->params.general.push_back({ (module.coords.node.data->begin() + i)->data(), 1, module.coords.node.isconst[i] ? 0 : info::mesh->dimension });
						insertNodeCoords = true;
					}
				}
				if (info::mesh->dimension == 3) {
					if (StringCompare::caseInsensitiveEq("COORDINATE_Z", value.evaluator[i * value.dimension + d]->variables[p])) {
						if (parameter.size.gp) {
							value.evaluator[i * value.dimension + d]->params.general.push_back({ (module.coords.gp.data->begin() + i)->data(), 2, module.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
							insertGpCoords = true;
						} else {
							value.evaluator[i * value.dimension + d]->params.general.push_back({ (module.coords.node.data->begin() + i)->data(), 2, module.coords.node.isconst[i] ? 0 : info::mesh->dimension });
							insertNodeCoords = true;
						}
					}
				}
				if (insertGpCoords) {
					parameter.addInput(module.coords.gp);
				}
				if (insertNodeCoords) {
					parameter.addInput(module.coords.node);
				}
//					if (StringCompare::caseInsensitiveEq("TIME", evaluators[i * dimension + d]->parameters[p])) {
//						info.addInput(_time, interval, dimension);
//					}
				if (value.evaluator[i * value.dimension + d]->params.general.size() == p) {
					eslog::error("ESPRESO internal error: implement dependency on parameter: '%s'\n", value.evaluator[i * value.dimension + d]->variables[p]);
				}
			}
		}
	}

	parameter.resize();

	module.addParameter(parameter);
	for (size_t i = 0; i < parameter.isconst.size(); ++i) {
		for (int d = 0; d < value.dimension && value.evaluator[i * value.dimension + d]; ++d) {
			module.elementOps[i].emplace_back(instantiate<AX_HeatTransfer::NGP, ExpressionsToParameter>(i, parameter, value.evaluator[i * value.dimension + d], d, value.dimension));
		}
	}
}

void fromExpression(AX_HeatTransfer &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &value)
{

}

}
