
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_EXPRESSION_H_

#include "basis/evaluator/evaluator.h"
#include "esinfo/meshinfo.h"
#include "mesh/mesh.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

class HeatTransferModuleOpt;

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

	template<int mask>
	ExpressionsToElements(ElementExternalParameter<mask> &parameter, double defaultValue, const char* name)
	: ExpressionsToParameter(parameter, defaultValue),
	  ElementOperatorBuilder(name),
	  parameter(parameter),
	  isset(parameter.isset)
	{
		parameter.builder = this;
	}

	bool build(HeatTransferModuleOpt &kernel) override;
	void apply(int interval);
};

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

	void replace(const std::string &name, BoundaryParameterPack &pack);

	bool build(HeatTransferModuleOpt &kernel) override;
	void apply(int region, int interval);
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
