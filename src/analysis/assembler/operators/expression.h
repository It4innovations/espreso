
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_

#include "basis/evaluator/evaluator.h"
#include "basis/utilities/parser.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "mesh/mesh.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"

#include <algorithm>

namespace espreso {

struct ExpressionsToParameter {

	int dimension;
	double defaultValue;
	std::vector<const Evaluator*> evaluators;
	std::vector<Evaluator::Params> params;

	template<int mask>
	ExpressionsToParameter(ElementExternalParameter<mask> &parameter, double defaultValue)
	: dimension(parameter.size.n * std::pow(info::mesh->dimension, parameter.size.ndimension)),
	  defaultValue(defaultValue),
	  evaluators(dimension * parameter.isconst.size(), NULL),
	  params(dimension * parameter.isconst.size())
	{

	}

	ExpressionsToParameter(BoundaryParameterData &parameter)
	: dimension(parameter.size.n * std::pow(info::mesh->dimension, parameter.size.ndimension)),
	  defaultValue(0), // dummy since boundary condition is always presented or there is no boundary condition
	  evaluators(dimension * parameter.isconst.size(), NULL),
	  params(dimension * parameter.isconst.size())
	{

	}
};

struct ExpressionsToElements: public ExpressionsToParameter, public ElementOperatorBuilder {
	ParameterData &parameter;
	std::vector<int> &isset;

	virtual ~ExpressionsToElements() {}

	template<int mask>
	ExpressionsToElements(ElementExternalParameter<mask> &parameter, double defaultValue, const char* name)
	: ExpressionsToParameter(parameter, defaultValue),
	  ElementOperatorBuilder(name),
	  parameter(parameter),
	  isset(parameter.isset)
	{
		parameter.builder = this;
	}
};

template<class Assembler>
struct ExpressionsToElementsSimple: public ExpressionsToElements {
	Assembler &assembler;

	template<int mask>
	ExpressionsToElementsSimple(Assembler &assembler, ElementExternalParameter<mask> &parameter, double defaultValue, const char* name)
	: ExpressionsToElements(parameter, defaultValue, name), assembler(assembler) {}

	bool build() override
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
							params[i * dimension + d].general.push_back({ (assembler.temp.initial.gp.data->begin() + i)->data(), 0, assembler.temp.initial.gp.isconst[i] ? 0 : 1 });
							parameter.addInput(assembler.temp.initial.gp);
						} else {
							params[i * dimension + d].general.push_back({ (assembler.temp.initial.node.data->begin() + i)->data(), 0, assembler.temp.initial.node.isconst[i] ? 0 : 1 });
							parameter.addInput(assembler.temp.initial.node);
						}
					}
					if (StringCompare::caseInsensitiveEq("TEMPERATURE", evaluators[i * dimension + d]->variables[p])) {
						if (parameter.size.gp) {
							params[i * dimension + d].general.push_back({ (assembler.temp.gp.data->begin() + i)->data(), 0, assembler.temp.gp.isconst[i] ? 0 : 1 });
							parameter.addInput(assembler.temp.gp);
						} else {
							params[i * dimension + d].general.push_back({ (assembler.temp.node.data->begin() + i)->data(), 0, assembler.temp.node.isconst[i] ? 0 : 1 });
							parameter.addInput(assembler.temp.node);
						}
					}
					bool insertGpCoords = false, insertNodeCoords = false;
					if (StringCompare::caseInsensitiveEq("COORDINATE_X", evaluators[i * dimension + d]->variables[p])) {
						if (parameter.size.gp) {
							params[i * dimension + d].general.push_back({ (assembler.coords.gp.data->begin() + i)->data(), 0, assembler.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
							insertGpCoords = true;
						} else {
							params[i * dimension + d].general.push_back({ (assembler.coords.node.data->begin() + i)->data(), 0, assembler.coords.node.isconst[i] ? 0 : info::mesh->dimension });
							insertNodeCoords = true;
						}
					}
					if (StringCompare::caseInsensitiveEq("COORDINATE_Y", evaluators[i * dimension + d]->variables[p])) {
						if (parameter.size.gp) {
							params[i * dimension + d].general.push_back({ (assembler.coords.gp.data->begin() + i)->data(), 1, assembler.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
							insertGpCoords = true;
						} else {
							params[i * dimension + d].general.push_back({ (assembler.coords.node.data->begin() + i)->data(), 1, assembler.coords.node.isconst[i] ? 0 : info::mesh->dimension });
							insertNodeCoords = true;
						}
					}
					if (info::mesh->dimension == 3) {
						if (StringCompare::caseInsensitiveEq("COORDINATE_Z", evaluators[i * dimension + d]->variables[p])) {
							if (parameter.size.gp) {
								params[i * dimension + d].general.push_back({ (assembler.coords.gp.data->begin() + i)->data(), 2, assembler.coords.gp.isconst[i] ? 0 : info::mesh->dimension });
								insertGpCoords = true;
							} else {
								params[i * dimension + d].general.push_back({ (assembler.coords.node.data->begin() + i)->data(), 2, assembler.coords.node.isconst[i] ? 0 : info::mesh->dimension });
								insertNodeCoords = true;
							}
						}
					}
					if (insertGpCoords) {
						parameter.addInput(assembler.coords.gp);
					}
					if (insertNodeCoords) {
						parameter.addInput(assembler.coords.node);
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

		assembler.addParameter(parameter);
		return true;
	}
	void apply(int interval)
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
};


//struct ExpressionsToElementsSimd: public ExpressionsToElements {
//
//	template<int mask>
//	ExpressionsToElementsSimd(ElementExternalParameter<mask> &parameter, double defaultValue, const char* name)
//	: ExpressionsToElements(parameter, defaultValue, name) {}
//
//	bool build(HeatTransferModuleOpt &assembler) override;
//	void apply(int interval);
//};

struct ExpressionsToBoundary: public BoundaryOperatorBuilder {
	BoundaryParameterPack &parameter;
	std::vector<ExpressionsToParameter> expressions;

	template<int mask>
	ExpressionsToBoundary(BoundaryExternalParameter<mask> &parameter, const char *name)
	: BoundaryOperatorBuilder(name), parameter(parameter)
	{
		expressions.reserve(info::mesh->boundaryRegions.size());
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			expressions.emplace_back(parameter.regions[r]);
		}
		parameter.builder = this;
	}

	void insert(int region, int dimension, const Evaluator *evaluator)
	{
		BoundaryParameterData &rdata = parameter.regions[region];
		rdata.isset = true;
		for (size_t i = 0; i < rdata.isconst.size(); ++i) {
			expressions[region].evaluators[i * expressions[region].dimension + dimension] = evaluator;
		}
	}

	void replace(const std::string &name, BoundaryParameterPack &pack)
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

	bool build() override
	{
		for (size_t r = 0; r < expressions.size(); ++r) {
//			ExpressionsToParameter &rexp = expressions[r];
//			BoundaryParameterData &rdata = parameter.regions[r];
//			for (size_t i = 0; i < rdata.isconst.size(); ++i) {
//				for (int d = 0; rexp.evaluators[i] && d < rexp.dimension; ++d) {
//					rdata.isset = true;
//					for (size_t p = 0; p < rexp.evaluators[i * rexp.dimension + d]->variables.size(); ++p) {
//						rdata.isconst[i] = false; // in the case of TIME it is possible to keep value constant
//						if (StringCompare::caseInsensitiveEq("INITIAL_TEMPERATURE", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
//							if (rdata.size.gp) {
//								rexp.params[i * rexp.dimension + d].general.push_back({ (assembler.temp.initial.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 0, assembler.temp.initial.boundary.gp.regions[rdata.region].isconst[i] ? 0 : 1 });
//								rdata.addInput(assembler.temp.initial.boundary.gp.regions[rdata.region]);
//							} else {
//								rexp.params[i * rexp.dimension + d].general.push_back({ (assembler.temp.initial.boundary.node.regions[rdata.region].data->begin() + i)->data(), 0, assembler.temp.initial.boundary.node.regions[rdata.region].isconst[i] ? 0 : 1 });
//								rdata.addInput(assembler.temp.initial.boundary.node.regions[rdata.region]);
//							}
//						}
//						if (StringCompare::caseInsensitiveEq("TEMPERATURE", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
//							if (rdata.size.gp) {
//								rexp.params[i * rexp.dimension + d].general.push_back({ (assembler.temp.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 0, assembler.temp.boundary.gp.regions[rdata.region].isconst[i] ? 0 : 1 });
//								rdata.addInput(assembler.temp.boundary.gp.regions[rdata.region]);
//							} else {
//								rexp.params[i * rexp.dimension + d].general.push_back({ (assembler.temp.boundary.node.regions[rdata.region].data->begin() + i)->data(), 0, assembler.temp.boundary.node.regions[rdata.region].isconst[i] ? 0 : 1 });
//								rdata.addInput(assembler.temp.boundary.node.regions[rdata.region]);
//							}
//						}
//						bool insertGpCoords = false, insertNodeCoords = false;
//						if (StringCompare::caseInsensitiveEq("COORDINATE_X", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
//							if (rdata.size.gp) {
//								rexp.params[i * rexp.dimension + d].general.push_back({ (assembler.coords.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 0, assembler.coords.boundary.gp.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
//								insertGpCoords = true;
//							} else {
//								rexp.params[i * rexp.dimension + d].general.push_back({ (assembler.coords.boundary.node.regions[rdata.region].data->begin() + i)->data(), 0, assembler.coords.boundary.node.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
//								insertNodeCoords = true;
//							}
//						}
//						if (StringCompare::caseInsensitiveEq("COORDINATE_Y", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
//							if (rdata.size.gp) {
//								rexp.params[i * rexp.dimension + d].general.push_back({ (assembler.coords.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 1, assembler.coords.boundary.gp.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
//								insertGpCoords = true;
//							} else {
//								rexp.params[i * rexp.dimension + d].general.push_back({ (assembler.coords.boundary.node.regions[rdata.region].data->begin() + i)->data(), 1, assembler.coords.boundary.node.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
//								insertNodeCoords = true;
//							}
//						}
//						if (info::mesh->dimension == 3) {
//							if (StringCompare::caseInsensitiveEq("COORDINATE_z", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
//								if (rdata.size.gp) {
//									rexp.params[i * rexp.dimension + d].general.push_back({ (assembler.coords.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 2, assembler.coords.boundary.gp.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
//									insertGpCoords = true;
//								} else {
//									rexp.params[i * rexp.dimension + d].general.push_back({ (assembler.coords.boundary.node.regions[rdata.region].data->begin() + i)->data(), 2, assembler.coords.boundary.node.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
//									insertNodeCoords = true;
//								}
//							}
//						}
//						if (insertGpCoords) {
//							rdata.addInput(assembler.coords.boundary.gp.regions[rdata.region]);
//						}
//						if (insertNodeCoords) {
//							rdata.addInput(assembler.coords.boundary.node.regions[rdata.region]);
//						}
//	//						if (StringCompare::caseInsensitiveEq("TIME", evaluators[i * dimension + d]->parameters[p])) {
//	//							info.addInput(_time, interval, dimension);
//	//						}
//						if (rexp.params[i * rexp.dimension + d].general.size() == p) {
//							eslog::error("ESPRESO internal error: implement dependency on parameter: '%s'\n", rexp.evaluators[i * expressions[r].dimension + d]->variables[p]);
//						}
//					}
//				}
//			}
//			rdata.resize();
//			assembler.addParameter(rdata);
		}
		return true;
	}
	void apply(int region, int interval)
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
};

struct ExpressionsToBoundaryFromElement: public ExpressionsToBoundary {
	template<int mask>
	ExpressionsToBoundaryFromElement(const ExpressionsToParameter &source, BoundaryExternalParameter<mask> &parameter, const char *name)
	: ExpressionsToBoundary(parameter, name)
	{
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				this->parameter.regions[r].isset = true;
			}
			for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
				for (int d = 0; d < expressions[r].dimension; ++d) {
					expressions[r].evaluators[i * expressions[r].dimension + d] = source.evaluators[info::mesh->boundaryRegions[r]->eintervals[i].einterval * expressions[r].dimension + d];
				}
			}
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_ */
