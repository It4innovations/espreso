
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_

#include "basis/evaluator/evaluator.h"
#include "basis/expression/variable.h"
#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"

namespace espreso {

template <class Setter>
struct ExternalExpression: ActionOperator {
	ExternalExpression(int interval, Evaluator *evaluator, const Setter &setter)
	: evaluator(evaluator),
	  params(evaluator->params),
	  setter(setter)
	{
		for (size_t i = 0; i < params.general.size(); ++i) {
			params.general[i].variable->set(interval, params.general[i]);
		}
		isconst = evaluator->variables.size() == 0;
	}

	Evaluator *evaluator;
	Evaluator::Params params;
	Setter setter;

	void move(int n)
	{
		for (size_t i = 0; i < params.general.size(); ++i) {
			params.general[i].val += n * params.general[i].increment;
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalNodeExpression: ExternalExpression<Setter>, Physics {
	using ExternalExpression<Setter>::ExternalExpression;

	void sisd(typename Physics::Element &element)
	{
		double results[nodes];
		this->evaluator->evalVector(nodes, this->params, results);
		for (size_t n = 0; n < nodes; ++n) {
			this->setter(element, n, 0, results[n]);
		}
		this->move(gps);
	}

	void simd(typename Physics::Element &element)
	{
		double results[SIMD::size * nodes];
		this->evaluator->evalVector(SIMD::size * nodes, this->params, results);
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				this->setter(element, n, s, results[nodes * s + n]); // TODO: check correct order
			}
		}
		this->move(SIMD::size * nodes);
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Setter>
struct ExternalGPsExpression: ExternalExpression<Setter>, Physics {
	using ExternalExpression<Setter>::ExternalExpression;

	void sisd(typename Physics::Element &element)
	{
		double results[gps];
		this->evaluator->evalVector(gps, this->params, results);
		for (size_t gp = 0; gp < gps; ++gp) {
			this->setter(element, gp, 0, results[gp]);
		}
		this->move(gps);
	}

	void simd(typename Physics::Element &element)
	{
		double results[SIMD::size * gps];
		this->evaluator->evalVector(SIMD::size * gps, this->params, results);
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				this->setter(element, gp, s, results[gps * s + gp]); // TODO: check correct order
			}
		}
		this->move(SIMD::size * gps);
	}
};

struct ExpressionsToParameter: public ActionOperator {
	ExpressionsToParameter(int interval, ParameterData &parameter, Evaluator *evaluator, size_t offset, size_t size)
	: data(parameter, interval),
	  evaluator(evaluator),
	  params(evaluator->params),
	  offset(offset), size(size)
	{
		for (size_t i = 0; i < params.general.size(); ++i) {
			params.general[i].variable->set(interval, params.general[i]);
		}
	}

	OutputParameterIterator data;
	Evaluator *evaluator;
	Evaluator::Params params;
	size_t offset, size;
};

template <size_t nodes, size_t gps>
struct ExpressionsToNodes: public ExpressionsToParameter {
	using ExpressionsToParameter::ExpressionsToParameter;

	void operator++()
	{
		++data;
		for (size_t i = 0; i < params.general.size(); ++i) {
			params.general[i].val += nodes * params.general[i].increment;
		}
	}

	void move(int n)
	{
		data += n;
		for (size_t i = 0; i < params.general.size(); ++i) {
			params.general[i].val += n * nodes * params.general[i].increment;
		}
	}

	void operator()()
	{
		double results[nodes];
		evaluator->evalVector(nodes, params, results);
		for (size_t n = 0; n < nodes; ++n) {
			data[n * size + offset] = results[n];
		}
	}
};

template <size_t nodes, size_t gps>
struct ExpressionsToGPs: public ExpressionsToParameter {
	using ExpressionsToParameter::ExpressionsToParameter;

	void operator++()
	{
		++data;
		for (size_t i = 0; i < params.general.size(); ++i) {
			params.general[i].val += gps * params.general[i].increment;
		}
	}

	void move(int n)
	{
		data += n;
		for (size_t i = 0; i < params.general.size(); ++i) {
			params.general[i].val += n * gps * params.general[i].increment;
		}
	}

	void operator()()
	{
		double results[gps];
		evaluator->evalVector(gps, params, results);
		for (size_t n = 0; n < gps; ++n) {
			data[n * size + offset] = results[n];
		}
	}
};

template<class Setter>
struct ExpressionsToParameter2: public ActionOperator {
	ExpressionsToParameter2(int interval, Setter setter, ParameterData &parameter, Evaluator *evaluator, size_t offset, size_t size)
	: setter(setter),
	  parameter(parameter, interval),
	  evaluator(evaluator),
	  params(evaluator->params),
	  offset(offset), size(size)
	{
		for (size_t i = 0; i < params.general.size(); ++i) {
			params.general[i].variable->set(interval, params.general[i]);
		}
	}

	Setter setter;
	OutputParameterIterator parameter;
	Evaluator *evaluator;
	Evaluator::Params params;
	size_t offset, size;
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, class Physics, class Setter>
struct ExpressionsToNodes2: public ExpressionsToParameter2<Setter>, Physics {
	using ExpressionsToParameter2<Setter>::ExpressionsToParameter;

	void operator++()
	{
		++this->parameter;
		for (size_t i = 0; i < this->params.general.size(); ++i) {
			this->params.general[i].val += nodes * this->params.general[i].increment;
		}
	}

	void move(int n)
	{
		this->parameter += n;
		for (size_t i = 0; i < this->params.general.size(); ++i) {
			this->params.general[i].val += n * nodes * this->params.general[i].increment;
		}
	}

	void operator()()
	{
		double results[nodes];
		this->evaluator->evalVector(nodes, this->params, results);
		for (size_t n = 0; n < nodes; ++n) {
//			this->parameter[n * this->size + this->offset] = results[n];
		}
	}

	void sisd(typename Physics::Element &element)
	{

	}

	void simd(typename Physics::Element &element)
	{
//		printf("expr simd to nodes\n");
		double results[SIMD::size * nodes];
		this->evaluator->evalVector(SIMD::size * nodes, this->params, results);
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t s = 0; s < SIMD::size; ++s) {
//				this->setter(element)[(n * this->size + this->offset) * SIMD::size + s] = results[s * SIMD::size + n];
			}
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, class Physics, class Setter>
struct ExpressionsToGPs2: public ExpressionsToParameter2<Setter>, Physics {
	using ExpressionsToParameter2<Setter>::ExpressionsToParameter2;

	void operator++()
	{
		++this->parameter;
		for (size_t i = 0; i < this->params.general.size(); ++i) {
			this->params.general[i].val += gps * this->params.general[i].increment;
		}
	}

	void move(int n)
	{
		this->parameter += n;
		for (size_t i = 0; i < this->params.general.size(); ++i) {
			this->params.general[i].val += n * gps * this->params.general[i].increment;
		}
	}

	void operator()()
	{
		double results[gps];
		this->evaluator->evalVector(gps, this->params, results);
		for (size_t n = 0; n < gps; ++n) {
//			this->parameter[n * this->size + this->offset] = results[n];
		}
	}

	void sisd(typename Physics::Element &element)
	{

	}

	void simd(typename Physics::Element &element)
	{
		double results[SIMD::size * gps];
		this->evaluator->evalVector(SIMD::size * gps, this->params, results);
		for (size_t n = 0; n < gps; ++n) {
			for (size_t s = 0; s < SIMD::size; ++s) {
//				this->setter(element)[(n * this->size + this->offset) * SIMD::size + s] = results[s * SIMD::size * n];
			}
		}
		move(SIMD::size);
	}
};

//struct ExpressionsToBoundary: public BoundaryOperatorBuilder {
//	BoundaryParameterPack &parameter;
//	std::vector<ExpressionsToParameter> expressions;
//
//	template<int mask>
//	ExpressionsToBoundary(BoundaryExternalParameter<mask> &parameter, const char *name)
//	: BoundaryOperatorBuilder(name), parameter(parameter)
//	{
////		expressions.reserve(info::mesh->boundaryRegions.size());
////		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
////			expressions.emplace_back(parameter.regions[r]);
////		}
////		parameter.builder = this;
//	}
//
//	void insert(int region, int dimension, const Evaluator *evaluator)
//	{
////		BoundaryParameterData &rdata = parameter.regions[region];
////		rdata.isset = true;
////		for (size_t i = 0; i < rdata.isconst.size(); ++i) {
////			expressions[region].evaluators[i * expressions[region].dimension + dimension] = evaluator;
////		}
//	}
//
//	void replace(const std::string &name, BoundaryParameterPack &pack)
//	{
////		for (size_t r = 0; r < expressions.size(); ++r) {
////			ExpressionsToParameter &rexp = expressions[r];
////			BoundaryParameterData &rdata = parameter.regions[r];
////			for (size_t i = 0; i < rdata.isconst.size(); ++i) {
////				for (int d = 0; rexp.evaluators[i] && d < rexp.dimension; ++d) {
////					for (size_t p = 0; p < rexp.evaluators[i * rexp.dimension + d]->variables.size(); ++p) {
////						if (StringCompare::caseInsensitiveEq(name, rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
////							rexp.params[i * rexp.dimension + d].general[p].val = (pack.regions[rdata.region].data->begin() + i)->data();
////							if (pack.regions[rdata.region].isconst[i]) {
////								rexp.params[i * rexp.dimension + d].general[p].increment = 0;
////								if (rexp.evaluators[i * rexp.dimension + d]->variables.size() == 1) {
////									rdata.isconst[i] = true;
////								}
////							}
////						}
////					}
////				}
////			}
////		}
//	}
//
//	void build()
//	{
//		for (size_t r = 0; r < expressions.size(); ++r) {
////			ExpressionsToParameter &rexp = expressions[r];
////			BoundaryParameterData &rdata = parameter.regions[r];
////			for (size_t i = 0; i < rdata.isconst.size(); ++i) {
////				for (int d = 0; rexp.evaluators[i] && d < rexp.dimension; ++d) {
////					rdata.isset = true;
////					for (size_t p = 0; p < rexp.evaluators[i * rexp.dimension + d]->variables.size(); ++p) {
////						rdata.isconst[i] = false; // in the case of TIME it is possible to keep value constant
////						if (StringCompare::caseInsensitiveEq("INITIAL_TEMPERATURE", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
////							if (rdata.size.gp) {
////								rexp.params[i * rexp.dimension + d].general.push_back({ (module.temp.initial.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 0, module.temp.initial.boundary.gp.regions[rdata.region].isconst[i] ? 0 : 1 });
////								rdata.addInput(module.temp.initial.boundary.gp.regions[rdata.region]);
////							} else {
////								rexp.params[i * rexp.dimension + d].general.push_back({ (module.temp.initial.boundary.node.regions[rdata.region].data->begin() + i)->data(), 0, module.temp.initial.boundary.node.regions[rdata.region].isconst[i] ? 0 : 1 });
////								rdata.addInput(module.temp.initial.boundary.node.regions[rdata.region]);
////							}
////						}
////						if (StringCompare::caseInsensitiveEq("TEMPERATURE", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
////							if (rdata.size.gp) {
////								rexp.params[i * rexp.dimension + d].general.push_back({ (module.temp.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 0, module.temp.boundary.gp.regions[rdata.region].isconst[i] ? 0 : 1 });
////								rdata.addInput(module.temp.boundary.gp.regions[rdata.region]);
////							} else {
////								rexp.params[i * rexp.dimension + d].general.push_back({ (module.temp.boundary.node.regions[rdata.region].data->begin() + i)->data(), 0, module.temp.boundary.node.regions[rdata.region].isconst[i] ? 0 : 1 });
////								rdata.addInput(module.temp.boundary.node.regions[rdata.region]);
////							}
////						}
////						bool insertGpCoords = false, insertNodeCoords = false;
////						if (StringCompare::caseInsensitiveEq("COORDINATE_X", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
////							if (rdata.size.gp) {
////								rexp.params[i * rexp.dimension + d].general.push_back({ (module.coords.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 0, module.coords.boundary.gp.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
////								insertGpCoords = true;
////							} else {
////								rexp.params[i * rexp.dimension + d].general.push_back({ (module.coords.boundary.node.regions[rdata.region].data->begin() + i)->data(), 0, module.coords.boundary.node.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
////								insertNodeCoords = true;
////							}
////						}
////						if (StringCompare::caseInsensitiveEq("COORDINATE_Y", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
////							if (rdata.size.gp) {
////								rexp.params[i * rexp.dimension + d].general.push_back({ (module.coords.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 1, module.coords.boundary.gp.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
////								insertGpCoords = true;
////							} else {
////								rexp.params[i * rexp.dimension + d].general.push_back({ (module.coords.boundary.node.regions[rdata.region].data->begin() + i)->data(), 1, module.coords.boundary.node.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
////								insertNodeCoords = true;
////							}
////						}
////						if (info::mesh->dimension == 3) {
////							if (StringCompare::caseInsensitiveEq("COORDINATE_z", rexp.evaluators[i * rexp.dimension + d]->variables[p])) {
////								if (rdata.size.gp) {
////									rexp.params[i * rexp.dimension + d].general.push_back({ (module.coords.boundary.gp.regions[rdata.region].data->begin() + i)->data(), 2, module.coords.boundary.gp.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
////									insertGpCoords = true;
////								} else {
////									rexp.params[i * rexp.dimension + d].general.push_back({ (module.coords.boundary.node.regions[rdata.region].data->begin() + i)->data(), 2, module.coords.boundary.node.regions[rdata.region].isconst[i] ? 0 : info::mesh->dimension });
////									insertNodeCoords = true;
////								}
////							}
////						}
////						if (insertGpCoords) {
////							rdata.addInput(module.coords.boundary.gp.regions[rdata.region]);
////						}
////						if (insertNodeCoords) {
////							rdata.addInput(module.coords.boundary.node.regions[rdata.region]);
////						}
////	//						if (StringCompare::caseInsensitiveEq("TIME", evaluators[i * dimension + d]->parameters[p])) {
////	//							info.addInput(_time, interval, dimension);
////	//						}
////						if (rexp.params[i * rexp.dimension + d].general.size() == p) {
////							eslog::error("ESPRESO internal error: implement dependency on parameter: '%s'\n", rexp.evaluators[i * expressions[r].dimension + d]->variables[p]);
////						}
////					}
////				}
////			}
////			rdata.resize();
////			module.addParameter(rdata);
//		}
//	}
//	void apply(int region, int interval)
//	{
////		if (parameter.regions[region].isset) {
////			if (parameter.regions[region].update[interval]) {
////				auto idata = parameter.regions[region].data->begin() + interval;
////				for (int d = 0; d < expressions[region].dimension; ++d) {
////					if (expressions[region].evaluators[expressions[region].dimension * interval + d]) {
////						expressions[region].evaluators[expressions[region].dimension * interval + d]->evalVector(idata->size() / expressions[region].dimension, expressions[region].dimension, expressions[region].params[expressions[region].dimension * interval + d], idata->data() + d);
////					}
////				}
////			}
////		}
//	}
//};
//
//struct ExpressionsToBoundaryFromElement: public ExpressionsToBoundary {
//	template<int mask>
//	ExpressionsToBoundaryFromElement(const ExpressionsToParameter &source, BoundaryExternalParameter<mask> &parameter, const char *name)
//	: ExpressionsToBoundary(parameter, name)
//	{
////		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
////			if (info::mesh->boundaryRegions[r]->dimension) {
////				this->parameter.regions[r].isset = true;
////			}
////			for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
////				for (int d = 0; d < expressions[r].dimension; ++d) {
////					expressions[r].evaluators[i * expressions[r].dimension + d] = source.evaluators[info::mesh->boundaryRegions[r]->eintervals[i].einterval * expressions[r].dimension + d];
////				}
////			}
////		}
//	}
//};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_ */
