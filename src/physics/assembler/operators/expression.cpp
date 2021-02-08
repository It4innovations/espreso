
#include "expression.h"

#include "physics/assembler/modules/heattransfer.module.opt.h"

#include "basis/utilities/parser.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

bool ExpressionsToElements::build(HeatTransferModuleOpt &kernel)
{
	if (ecfname.size() == 0) {
		return false;
	}
	for (size_t i = 0; i < parameter.isconst.size(); ++i) {
		for (int d = 0; d < dimension && evaluators[i * dimension + d]; ++d) {
			for (size_t p = 0; p < evaluators[i * dimension + d]->variables.size(); ++p) {
				parameter.isconst[i] = false; // in the case of TIME it is possible to keep value constant
				if (StringCompare::caseInsensitiveEq("INITIAL_TEMPERATURE", evaluators[i * dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						params[i * dimension + d].general.push_back({ (kernel.temp.initial.gp.data->begin() + i)->data(), 0, kernel.temp.initial.gp.isconst[i] ? 0 : 1 });
					} else {
						params[i * dimension + d].general.push_back({ (kernel.temp.initial.node.data->begin() + i)->data(), 0, kernel.temp.initial.node.isconst[i] ? 0 : 1 });
					}
				}
				if (StringCompare::caseInsensitiveEq("TEMPERATURE", evaluators[i * dimension + d]->variables[p])) {
					update[i] += 1;
					if (parameter.size.gp) {
						params[i * dimension + d].general.push_back({ (kernel.temp.gp.data->begin() + i)->data(), 0, kernel.temp.gp.isconst[i] ? 0 : 1 });
					} else {
						params[i * dimension + d].general.push_back({ (kernel.temp.node.data->begin() + i)->data(), 0, kernel.temp.node.isconst[i] ? 0 : 1 });
					}
				}
				if (StringCompare::caseInsensitiveEq("COORDINATE_X", evaluators[i * dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						params[i * dimension + d].general.push_back({ (kernel.coords.gp.data->begin() + i)->data(), 0, kernel.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
					} else {
						params[i * dimension + d].general.push_back({ (kernel.coords.node.data->begin() + i)->data(), 0, kernel.coords.node.isconst[i] ? 0 : info::mesh->dimension });
					}
				}
				if (StringCompare::caseInsensitiveEq("COORDINATE_Y", evaluators[i * dimension + d]->variables[p])) {
					if (parameter.size.gp) {
						params[i * dimension + d].general.push_back({ (kernel.coords.gp.data->begin() + i)->data(), 1, kernel.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
					} else {
						params[i * dimension + d].general.push_back({ (kernel.coords.node.data->begin() + i)->data(), 1, kernel.coords.node.isconst[i] ? 0 : info::mesh->dimension });
					}
				}
				if (info::mesh->dimension == 3) {
					if (StringCompare::caseInsensitiveEq("COORDINATE_Z", evaluators[i * dimension + d]->variables[p])) {
						if (parameter.size.gp) {
							params[i * dimension + d].general.push_back({ (kernel.coords.gp.data->begin() + i)->data(), 2, kernel.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
						} else {
							params[i * dimension + d].general.push_back({ (kernel.coords.node.data->begin() + i)->data(), 2, kernel.coords.node.isconst[i] ? 0 : info::mesh->dimension });
						}
					}
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
	return true;
}

void ExpressionsToElements::apply(int interval)
{
	auto idata = parameter.data->begin() + interval;
	if (update[interval] || !parameter.version[interval]) {
		if (parameter.isconst[interval]) {
			if (Operator::print) printf("\tOP::ELEMENT::%d::EXPR::CONST\n", interval);
		} else {
			if (Operator::print) printf("\tOP::ELEMENT::%d::EXPR::EVAL\n", interval);
		}

		for (int d = 0; d < dimension; ++d) {
			if (evaluators[dimension * interval + d]) {
				evaluators[dimension * interval + d]->evalVector(idata->size() / dimension, dimension, params[dimension * interval + d], idata->data() + d);
			} else {
				for (size_t j = 0; j < idata->size(); j += dimension) {
					idata->data()[j + d] = defaultValue;
				}
			}
		}
		++parameter.version[interval];
	} else {
		if (Operator::print > 2) printf("\tOP::ELEMENT::%d::EXPR::SKIPPED\n", interval);
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
		rdata.resize();
	}
}

bool ExpressionsToBoundary::build(HeatTransferModuleOpt &kernel)
{
	if (ecfname.size() == 0) {
		return false;
	}
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
						} else {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.temp.initial.boundary.node.regions[rdata.region].data->begin() + i)->data(), 0, kernel.temp.initial.boundary.node.regions[rdata.region].isconst[i] ? 0 : 1 });
						}
					}
					if (StringCompare::caseInsensitiveEq("TEMPERATURE", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
						rexp.update[i] += 1;
						if (rdata.size.gp) {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.temp.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 0, kernel.temp.boundary.gp.regions[rdata.region].isconst[i] ? 0 : 1 });
						} else {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.temp.boundary.node.regions[rdata.region].data->begin() + i)->data(), 0, kernel.temp.boundary.node.regions[rdata.region].isconst[i] ? 0 : 1 });
						}
					}
					if (StringCompare::caseInsensitiveEq("COORDINATE_X", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
						if (rdata.size.gp) {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.coords.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 0, kernel.coords.boundary.gp.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
						} else {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.coords.boundary.node.regions[rdata.region].data->begin() + i)->data(), 0, kernel.coords.boundary.node.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
						}
					}
					if (StringCompare::caseInsensitiveEq("COORDINATE_Y", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
						if (rdata.size.gp) {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.coords.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 1, kernel.coords.boundary.gp.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
						} else {
							rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.coords.boundary.node.regions[rdata.region].data->begin() + i)->data(), 1, kernel.coords.boundary.node.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
						}
					}
					if (info::mesh->dimension == 3) {
						if (StringCompare::caseInsensitiveEq("COORDINATE_z", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
							if (rdata.size.gp) {
								rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.coords.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 2, kernel.coords.boundary.gp.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
							} else {
								rexp.params[i * rexp.dimension + d].general.push_back({ (kernel.coords.boundary.node.regions[rdata.region].data->begin() + i)->data(), 2, kernel.coords.boundary.node.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
							}
						}
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
	}
	return true;
}

void ExpressionsToBoundary::apply(int region, int interval)
{
	if (parameter.regions[region].isset) {
		if (expressions[region].update[interval] || !parameter.regions[region].version[interval]) {
			if (parameter.regions[region].isconst[interval]) {
				if (Operator::print) printf("\tOP::BOUNDARY::%d::%d::EXPR::CONST\n", region, interval);
			} else {
				if (Operator::print) printf("\tOP::BOUNDARY::%d::%d::EXPR::EVAL\n", region, interval);
			}

			auto idata = parameter.regions[region].data->begin() + interval;
			for (int d = 0; d < expressions[region].dimension; ++d) {
				if (expressions[region].evaluators[expressions[region].dimension * interval + d]) {
					expressions[region].evaluators[expressions[region].dimension * interval + d]->evalVector(idata->size() / expressions[region].dimension, expressions[region].dimension, expressions[region].params[expressions[region].dimension * interval + d], idata->data() + d);
				}
			}
			++parameter.regions[region].version[interval];
		} else {
			if (Operator::print > 2) printf("\tOP::BOUNDARY::%d::%d::EXPR::SKIPPED\n", region, interval);
		}
	}
}
