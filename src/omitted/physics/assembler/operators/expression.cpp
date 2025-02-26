
#include "expression.h"

#include "physics/assembler/modules/heattransfer.module.opt.h"

#include "basis/utilities/parser.h"
#include "esinfo/eslog.hpp"

#include <algorithm>

using namespace espreso;

bool ExpressionsToElementsSimple::build(HeatTransferModuleOpt &kernel)
{
	if (std::all_of(evaluators.begin(), evaluators.end(), [] (const Evaluator *ev) { return ev == NULL; })) {
		return false;
	}
	for (size_t i = 0; i < parameter.isconst.size(); ++i) {
		for (int d = 0; d < dimension && evaluators[i * dimension + d]; ++d) {
			isset[i] = isset[i] | evaluators[i * dimension + d]->isset;
			for (size_t p = 0; p < evaluators[i * dimension + d]->variables.size(); ++p) {
				parameter.isconst[i] = false; // in the case of TIME it is possible to keep value constant
				if (StringCompare::caseInsensitiveEq("INITIAL_TEMPERATURE", evaluators[i * dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						params[i * dimension + d].general.push_back({ (kernel.temp.initial.gp.data->begin() + i)->data(), 0, kernel.temp.initial.gp.isconst[i] ? 0 : 1 });
						parameter.addInput(kernel.temp.initial.gp);
					} else {
						params[i * dimension + d].general.push_back({ (kernel.temp.initial.node.data->begin() + i)->data(), 0, kernel.temp.initial.node.isconst[i] ? 0 : 1 });
						parameter.addInput(kernel.temp.initial.node);
					}
				}
				if (StringCompare::caseInsensitiveEq("TEMPERATURE", evaluators[i * dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						params[i * dimension + d].general.push_back({ (kernel.temp.gp.data->begin() + i)->data(), 0, kernel.temp.gp.isconst[i] ? 0 : 1 });
						parameter.addInput(kernel.temp.gp);
					} else {
						params[i * dimension + d].general.push_back({ (kernel.temp.node.data->begin() + i)->data(), 0, kernel.temp.node.isconst[i] ? 0 : 1 });
						parameter.addInput(kernel.temp.node);
					}
				}
				bool insertGpCoords = false, insertNodeCoords = false;
				if (StringCompare::caseInsensitiveEq("COORDINATE_X", evaluators[i * dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						params[i * dimension + d].general.push_back({ (kernel.coords.gp.data->begin() + i)->data(), 0, kernel.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
						insertGpCoords = true;
					} else {
						params[i * dimension + d].general.push_back({ (kernel.coords.node.data->begin() + i)->data(), 0, kernel.coords.node.isconst[i] ? 0 : info::mesh->dimension });
						insertNodeCoords = true;
					}
				}
				if (StringCompare::caseInsensitiveEq("COORDINATE_Y", evaluators[i * dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						params[i * dimension + d].general.push_back({ (kernel.coords.gp.data->begin() + i)->data(), 1, kernel.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
						insertGpCoords = true;
					} else {
						params[i * dimension + d].general.push_back({ (kernel.coords.node.data->begin() + i)->data(), 1, kernel.coords.node.isconst[i] ? 0 : info::mesh->dimension });
						insertNodeCoords = true;
					}
				}
				if (info::mesh->dimension == 3) {
					if (StringCompare::caseInsensitiveEq("COORDINATE_Z", evaluators[i * dimension + d]->variables[p])) {
						if (parameter.size.gp) {
							params[i * dimension + d].general.push_back({ (kernel.coords.gp.data->begin() + i)->data(), 2, kernel.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
							insertGpCoords = true;
						} else {
							params[i * dimension + d].general.push_back({ (kernel.coords.node.data->begin() + i)->data(), 2, kernel.coords.node.isconst[i] ? 0 : info::mesh->dimension });
							insertNodeCoords = true;
						}
					}
				}
				if (insertGpCoords) {
					parameter.addInput(kernel.coords.gp);
				}
				if (insertNodeCoords) {
					parameter.addInput(kernel.coords.node);
				}
//					if (StringCompare::caseInsensitiveEq("TIME", evaluators[i * dimension + d]->parameters[p])) {
//						info.addInput(_time, interval, dimension);
//					}
				if (params[i * dimension + d].general.size() == p) {
					eslog::error("ESPRESO internal error: implement dependency on parameter: '%s'\n", evaluators[i * dimension + d]->variables[p]);
				}
			}
		}
	}
	
	parameter.resize();
	
	kernel.addParameter(parameter);
	return true;
}

void ExpressionsToElementsSimple::apply(int interval)
{
	auto idata = parameter.data->begin() + interval;
	if (parameter.update[interval]) {
		for (int d = 0; d < dimension; ++d) {
			if (evaluators[dimension * interval + d]) {
				evaluators[dimension * interval + d]->evalVector(idata->size() / dimension, dimension, params[dimension * interval + d], idata->data() + d);
			} else {
				for (size_t j = 0; j < idata->size(); j += dimension) {
					idata->data()[j + d] = defaultValue;
				}
			}
		}
	}
}

bool ExpressionsToElementsSimd::build(HeatTransferModuleOpt &kernel)
{
	if (std::all_of(evaluators.begin(), evaluators.end(), [] (const Evaluator *ev) { return ev == NULL; })) {
		return false;
	}
	for (size_t i = 0; i < parameter.isconst.size(); ++i) {
		for (int d = 0; d < dimension && evaluators[i * dimension + d]; ++d) {
			isset[i] = isset[i] | evaluators[i * dimension + d]->isset;
			for (size_t p = 0; p < evaluators[i * dimension + d]->variables.size(); ++p) {
				parameter.isconst[i] = false; // in the case of TIME it is possible to keep value constant
				if (StringCompare::caseInsensitiveEq("INITIAL_TEMPERATURE", evaluators[i * dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						params[i * dimension + d].general.push_back({ (kernel.temp.initial.gp.data->begin() + i)->data(), 0, kernel.temp.initial.gp.isconst[i] ? 0 : 1 });
						parameter.addInput(kernel.temp.initial.gp);
					} else {
						params[i * dimension + d].general.push_back({ (kernel.temp.initial.node.data->begin() + i)->data(), 0, kernel.temp.initial.node.isconst[i] ? 0 : 1 });
						parameter.addInput(kernel.temp.initial.node);
					}
				}
				if (StringCompare::caseInsensitiveEq("TEMPERATURE", evaluators[i * dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						params[i * dimension + d].general.push_back({ (kernel.temp.gp.data->begin() + i)->data(), 0, kernel.temp.gp.isconst[i] ? 0 : 1 });
						parameter.addInput(kernel.temp.gp);
					} else {
						params[i * dimension + d].general.push_back({ (kernel.temp.node.data->begin() + i)->data(), 0, kernel.temp.node.isconst[i] ? 0 : 1 });
						parameter.addInput(kernel.temp.node);
					}
				}
				bool insertGpCoords = false, insertNodeCoords = false;
				if (StringCompare::caseInsensitiveEq("COORDINATE_X", evaluators[i * dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						params[i * dimension + d].general.push_back({ (kernel.coordsSimd.gp.data->begin() + i)->data(), 0, kernel.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
						insertGpCoords = true;
					} else {
						params[i * dimension + d].general.push_back({ (kernel.coordsSimd.node.data->begin() + i)->data(), 0, kernel.coords.node.isconst[i] ? 0 : info::mesh->dimension });
						insertNodeCoords = true;
					}
				}
				if (StringCompare::caseInsensitiveEq("COORDINATE_Y", evaluators[i * dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						params[i * dimension + d].general.push_back({ (kernel.coordsSimd.gp.data->begin() + i)->data(), 1, kernel.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
						insertGpCoords = true;
					} else {
						params[i * dimension + d].general.push_back({ (kernel.coordsSimd.node.data->begin() + i)->data(), 1, kernel.coords.node.isconst[i] ? 0 : info::mesh->dimension });
						insertNodeCoords = true;
					}
				}
				if (info::mesh->dimension == 3) {
					if (StringCompare::caseInsensitiveEq("COORDINATE_Z", evaluators[i * dimension + d]->variables[p])) {
						if (parameter.size.gp) {
							params[i * dimension + d].general.push_back({ (kernel.coordsSimd.gp.data->begin() + i)->data(), 2, kernel.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
							insertGpCoords = true;
						} else {
							params[i * dimension + d].general.push_back({ (kernel.coordsSimd.node.data->begin() + i)->data(), 2, kernel.coords.node.isconst[i] ? 0 : info::mesh->dimension });
							insertNodeCoords = true;
						}
					}
				}
				if (insertGpCoords) {
					parameter.addInput(kernel.coords.gp);
				}
				if (insertNodeCoords) {
					parameter.addInput(kernel.coords.node);
				}
//					if (StringCompare::caseInsensitiveEq("TIME", evaluators[i * dimension + d]->parameters[p])) {
//						info.addInput(_time, interval, dimension);
//					}
				if (params[i * dimension + d].general.size() == p) {
					eslog::error("ESPRESO internal error: implement dependency on parameter: '%s'\n", evaluators[i * dimension + d]->variables[p]);
				}
			}
		}
	}
	
	parameter.resizeAligned(SIMD::size*sizeof(double));
	
	kernel.addParameter(parameter);
	return true;
}

void ExpressionsToElementsSimd::apply(int interval)
{
	auto idata = parameter.data->begin() + interval;
	if (parameter.update[interval]) {
		for (int d = 0; d < dimension; ++d) {
			if (evaluators[dimension * interval + d]) {
				evaluators[dimension * interval + d]->evalVectorSimd(idata->size() / dimension, dimension, params[dimension * interval + d], idata->data() + d);
			} else {
				for (size_t j = 0; j < idata->size(); j += dimension) {
					idata->data()[j + d] = defaultValue;
				}
			}
		}
	}
}

void ExpressionsToBoundary::replace(const std::string &name, BoundaryParameterPack &pack)
{
	for (size_t r = 0; r < expressions.size(); ++r) {
		ExpressionsToParameter &rexp = expressions[r];
		BoundaryParameterData &rdata = parameter.regions[r];
		for (size_t i = 0; i < rdata.isconst.size(); ++i) {
			for (int d = 0; rexp.evaluators[i] && d < rexp.dimension; ++d) {
				for (size_t p = 0; p < rexp.evaluators[i * rexp.dimension + d]->variables.size(); ++p) {
					if (StringCompare::caseInsensitiveEq(name, rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
						rexp.params[i * rexp.dimension + d].general[p].val = (pack.regions[rdata.region].data->begin() + i)->data();
						if (pack.regions[rdata.region].isconst[i]) {
							rexp.params[i * rexp.dimension + d].general[p].increment = 0;
							if (rexp.evaluators[i * rexp.dimension + d]->variables.size() == 1) {
								rdata.isconst[i] = true;
							}
						}
					}
				}
			}
		}
	}
}

bool ExpressionsToBoundary::build(HeatTransferModuleOpt &kernel)
{
	for (size_t r = 0; r < expressions.size(); ++r) {
		ExpressionsToParameter &rexp = expressions[r];
		BoundaryParameterData &rdata = parameter.regions[r];
		for (size_t i = 0; i < rdata.isconst.size(); ++i) {
			for (int d = 0; rexp.evaluators[i] && d < rexp.dimension; ++d) {
				rdata.isset = true;
				for (size_t p = 0; p < rexp.evaluators[i * rexp.dimension + d]->variables.size(); ++p) {
					rdata.isconst[i] = false; // in the case of TIME it is possible to keep value constant
					if (StringCompare::caseInsensitiveEq("INITIAL_TEMPERATURE", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
						if (rdata.size.gp) {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.temp.initial.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 0, kernel.temp.initial.boundary.gp.regions[rdata.region].isconst[i] ? 0 : 1 });
							rdata.addInput(kernel.temp.initial.boundary.gp.regions[rdata.region]);
						} else {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.temp.initial.boundary.node.regions[rdata.region].data->begin() + i)->data(), 0, kernel.temp.initial.boundary.node.regions[rdata.region].isconst[i] ? 0 : 1 });
							rdata.addInput(kernel.temp.initial.boundary.node.regions[rdata.region]);
						}
					}
					if (StringCompare::caseInsensitiveEq("TEMPERATURE", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
						if (rdata.size.gp) {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.temp.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 0, kernel.temp.boundary.gp.regions[rdata.region].isconst[i] ? 0 : 1 });
							rdata.addInput(kernel.temp.boundary.gp.regions[rdata.region]);
						} else {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.temp.boundary.node.regions[rdata.region].data->begin() + i)->data(), 0, kernel.temp.boundary.node.regions[rdata.region].isconst[i] ? 0 : 1 });
							rdata.addInput(kernel.temp.boundary.node.regions[rdata.region]);
						}
					}
					bool insertGpCoords = false, insertNodeCoords = false;
					if (StringCompare::caseInsensitiveEq("COORDINATE_X", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
						if (rdata.size.gp) {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.coords.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 0, kernel.coords.boundary.gp.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
							insertGpCoords = true;
						} else {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.coords.boundary.node.regions[rdata.region].data->begin() + i)->data(), 0, kernel.coords.boundary.node.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
							insertNodeCoords = true;
						}
					}
					if (StringCompare::caseInsensitiveEq("COORDINATE_Y", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
						if (rdata.size.gp) {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.coords.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 1, kernel.coords.boundary.gp.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
							insertGpCoords = true;
						} else {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.coords.boundary.node.regions[rdata.region].data->begin() + i)->data(), 1, kernel.coords.boundary.node.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
							insertNodeCoords = true;
						}
					}
					if (info::mesh->dimension == 3) {
						if (StringCompare::caseInsensitiveEq("COORDINATE_z", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
							if (rdata.size.gp) {
								rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.coords.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 2, kernel.coords.boundary.gp.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
								insertGpCoords = true;
							} else {
								rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.coords.boundary.node.regions[rdata.region].data->begin() + i)->data(), 2, kernel.coords.boundary.node.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
								insertNodeCoords = true;
							}
						}
					}
					if (insertGpCoords) {
						rdata.addInput(kernel.coords.boundary.gp.regions[rdata.region]);
					}
					if (insertNodeCoords) {
						rdata.addInput(kernel.coords.boundary.node.regions[rdata.region]);
					}
//						if (StringCompare::caseInsensitiveEq("TIME", evaluators[i * dimension + d]->parameters[p])) {
//							info.addInput(_time, interval, dimension);
//						}
					if (rexp.params[i * rexp.dimension + d].general.size() == p) {
						eslog::error("ESPRESO internal error: implement dependency on parameter: '%s'\n", rexp.evaluators[i * expressions[r].dimension + d]->variables[p]);
					}
				}
			}
		}
		rdata.resize();
		kernel.addParameter(rdata);
	}
	return true;
}

void ExpressionsToBoundary::apply(int region, int interval)
{
	if (parameter.regions[region].isset) {
		if (parameter.regions[region].update[interval]) {
			auto idata = parameter.regions[region].data->begin() + interval;
			for (int d = 0; d < expressions[region].dimension; ++d) {
				if (expressions[region].evaluators[expressions[region].dimension * interval + d]) {
					expressions[region].evaluators[expressions[region].dimension * interval + d]->evalVector(idata->size() / expressions[region].dimension, expressions[region].dimension, expressions[region].params[expressions[region].dimension * interval + d], idata->data() + d);
				}
			}
		}
	}
}
