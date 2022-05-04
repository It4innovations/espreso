
#include "variable.h"
#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "basis/evaluator/expressionevaluator.h"
#include "basis/expression/expression.h"
#include "basis/utilities/parser.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "config/holders/expression.h"
#include "mesh/store/nameddata.h"

decltype(espreso::Variable::list) espreso::Variable::list;

using namespace espreso;

void Variable::gather(size_t regions)
{
	for (auto range = info::ecf->ranges.begin(); range != info::ecf->ranges.end(); ++range) {
		Variable::list.global[range->first] = new RangeVariable(range->second);
	}

	for (auto it = ECFExpression::parametrized.begin(); it != ECFExpression::parametrized.end(); ++it) {
		std::string &expr = (*it)->value;
		for (size_t i = 0; i + 1 < expr.size(); ++i) {
			if (expr[i] == ':' && expr[i + 1] == ':') {
				expr[i] = expr[i + 1] = '_';
			}
		}

		std::vector<std::string> variables;
		Expression::collectVariables(expr, variables);

		for (size_t i = 0; i < variables.size(); ++i) {
			std::string variable = Parser::uppercase(variables[i]);

			if (Variable::list.global.find(variable) != Variable::list.global.end()) {
				(*it)->parameters.push_back(variable);
				continue;
			}

			if (StringCompare::caseInsensitivePreffix("GLOBAL__", variable)) {
				Variable::list.global.insert(std::make_pair(variable.substr(8), nullptr));
				(*it)->parameters.push_back(variable.substr(8));
				continue;
			}
			if (StringCompare::caseInsensitivePreffix("ELEMENT__", variable)) {
				Variable::list.element.insert(std::make_pair(variable.substr(9), nullptr));
				(*it)->parameters.push_back(variable.substr(9));
				continue;
			}
			if (StringCompare::caseInsensitivePreffix("ENODES__", variable)) {
				Variable::list.enodes.insert(std::make_pair(variable.substr(8), nullptr));
				(*it)->parameters.push_back(variable.substr(8));
				continue;
			}
			if (StringCompare::caseInsensitivePreffix("EGPS__", variable)) {
				Variable::list.egps.insert(std::make_pair(variable.substr(6), nullptr));
				(*it)->parameters.push_back(variable.substr(6));
				continue;
			}
//			if (StringCompare::caseInsensitivePreffix("BNODES__", variable)) {
//				Variable::list.region.front().enodes.insert(std::make_pair(variable.substr(8), nullptr));
//				(*it)->parameters.push_back(variable.substr(8));
//				continue;
//			}
//			if (StringCompare::caseInsensitivePreffix("BGPS__", variable)) {
//				Variable::list.region.front().egps.insert(std::make_pair(variable.substr(6), nullptr));
//				(*it)->parameters.push_back(variable.substr(6));
//				continue;
//			}
			if (StringCompare::caseInsensitivePreffix("NODE__", variable)) {
				Variable::list.node.insert(std::make_pair(variable.substr(6), nullptr));
				(*it)->parameters.push_back(variable.substr(6));
				continue;
			}
			switch ((*it)->scope) {
			case ECFExpression::Scope::GLOBAL: Variable::list.global.insert(std::make_pair(variable, nullptr)); break;
			case ECFExpression::Scope::ELEMENT: Variable::list.element.insert(std::make_pair(variable, nullptr)); break;
			case ECFExpression::Scope::ENODES: Variable::list.enodes.insert(std::make_pair(variable, nullptr)); break;
			case ECFExpression::Scope::EGPS: Variable::list.egps.insert(std::make_pair(variable, nullptr)); break;
//			case ECFExpression::Scope::BNODES: Variable::list.region.front().enodes.insert(std::make_pair(variable, nullptr)); break;
//			case ECFExpression::Scope::BGPS: Variable::list.region.front().egps.insert(std::make_pair(variable, nullptr)); break;
			case ECFExpression::Scope::NODE: Variable::list.node.insert(std::make_pair(variable, nullptr)); break;
			default: continue;
			}
			(*it)->parameters.push_back(variable);
		}
	}
	Variable::list.region.resize(regions);
}

void Variable::analyze(ECFExpression &expr, size_t region)
{
	std::vector<std::string> variables;
	Expression::collectVariables(expr.value, variables);

	for (size_t i = 0; i < variables.size(); ++i) {
		std::string variable = Parser::uppercase(variables[i]);
		if (StringCompare::caseInsensitivePreffix("BNODES__", variable)) {
			Variable::list.region[region].enodes.insert(std::make_pair(variable.substr(8), nullptr));
			expr.parameters.push_back(variable.substr(8));
			continue;
		}
		if (StringCompare::caseInsensitivePreffix("BGPS__", variable)) {
			Variable::list.region[region].egps.insert(std::make_pair(variable.substr(6), nullptr));
			expr.parameters.push_back(variable.substr(6));
			continue;
		}
		if (!Variable::list.global.count(variable)) {
			switch (expr.scope) {
			case ECFExpression::Scope::BNODES: Variable::list.region[region].enodes.insert(std::make_pair(variable, nullptr)); break;
			case ECFExpression::Scope::BGPS: Variable::list.region[region].egps.insert(std::make_pair(variable, nullptr)); break;
			default: continue;
			}
			expr.parameters.push_back(variable);
		}
	}
}


bool Variable::create(ECFExpression &expr)
{
	std::vector<std::string> variables;
	Expression::collectVariables(expr.value, variables);

	Evaluator::Params params;

	for (size_t i = 0; i < variables.size(); ++i) {
		std::map<std::string, Variable*> *map = nullptr;
		std::string variable = Parser::uppercase(variables[i]);

		if (StringCompare::caseInsensitivePreffix("GLOBAL__", variable)) {
			map = &Variable::list.global;
			variable = variable.substr(8);
		}
		if (StringCompare::caseInsensitivePreffix("ELEMENT__", variable)) {
			map = &Variable::list.element;
			variable = variable.substr(9);
		}
		if (StringCompare::caseInsensitivePreffix("ENODES__", variable)) {
			map = &Variable::list.enodes;
			variable = variable.substr(8);
		}
		if (StringCompare::caseInsensitivePreffix("EGPS__", variable)) {
			map = &Variable::list.egps;
			variable = variable.substr(6);
		}
		if (StringCompare::caseInsensitivePreffix("NODE__", variable)) {
			map = &Variable::list.node;
			variable = variable.substr(6);
		}

		if (Variable::list.global.find(variable) != Variable::list.global.end()) {
			params.general.push_back(Evaluator::Params::General{nullptr, 0, 0, Variable::list.global.find(variable)->second});
			continue;
		}

		switch (expr.scope) {
		case ECFExpression::Scope::GLOBAL: map = &Variable::list.global; break;
		case ECFExpression::Scope::ELEMENT: map = &Variable::list.element; break;
		case ECFExpression::Scope::ENODES: map = &Variable::list.enodes; break;
		case ECFExpression::Scope::EGPS: map = &Variable::list.egps; break;
		case ECFExpression::Scope::NODE: map = &Variable::list.node; break;
		default: break;
		}

		std::map<std::string, Variable*>::iterator it;
		if (map == nullptr) {
			if ((it = Variable::list.global.find(variable)) == Variable::list.global.end()) {
				return false;
			}
		} else {
			if ((it = map->find(variable)) == map->end()) {
				if ((it = Variable::list.global.find(variable)) == Variable::list.global.end()) {
					return false;
				}
			}
		}
		params.general.push_back(Evaluator::Params::General{nullptr, 0, 0, it->second}); // we will set parameter later for each interval
	}
	expr.evaluator = new ExpressionEvaluator(expr.value, expr.parameters);
	expr.evaluator->params = params;

	return true;
}

bool Variable::create(ECFExpression &expr, size_t region)
{
	std::vector<std::string> variables;
	Expression::collectVariables(expr.value, variables);

	Evaluator::Params params;

	for (size_t i = 0; i < variables.size(); ++i) {
		std::map<std::string, Variable*> *map = nullptr;
		std::string variable = Parser::uppercase(variables[i]);

		if (StringCompare::caseInsensitivePreffix("GLOBAL__", variable)) {
			map = &Variable::list.global;
			variable = variable.substr(8);
		}
		if (StringCompare::caseInsensitivePreffix("BNODES__", variable)) {
			map = &Variable::list.region[region].enodes;
			variable = variable.substr(8);
		}
		if (StringCompare::caseInsensitivePreffix("BGPS__", variable)) {
			map = &Variable::list.region[region].egps;
			variable = variable.substr(6);
		}
		if (StringCompare::caseInsensitivePreffix("NODE__", variable)) {
			map = &Variable::list.node;
			variable = variable.substr(6);
		}
		switch (expr.scope) {
		case ECFExpression::Scope::GLOBAL: map = &Variable::list.global; break;
		case ECFExpression::Scope::BNODES: map = &Variable::list.region[region].enodes; break;
		case ECFExpression::Scope::BGPS: map = &Variable::list.region[region].egps; break;
		case ECFExpression::Scope::NODE: map = &Variable::list.node; break;
		default: break;
		}

		std::map<std::string, Variable*>::iterator it;
		if (map == nullptr) {
			if ((it = Variable::list.global.find(variable)) == Variable::list.global.end()) {
				return false;
			}
		} else {
			if ((it = map->find(variable)) == map->end()) {
				if ((it = Variable::list.global.find(variable)) == Variable::list.global.end()) {
					return false;
				}
			}
		}
		params.general.push_back(Evaluator::Params::General{nullptr, 0, 0, it->second});
	}
	expr.evaluator = new ExpressionEvaluator(expr.value, expr.parameters);
	expr.evaluator->params = params;

	return true;
}

void Variable::print()
{
	auto _print = [] (std::map<std::string, Variable*> &map) {
		for (auto it = map.begin(); it != map.end(); ++it) {
			eslog::info(" =   %s %*s = \n", it->first.c_str(), 86 - it->first.size(), "");
		}
	};

	eslog::info(" = GLOBAL VARIABLES ========================================================================== \n");
	_print(list.global);
	eslog::info(" = PER ELEMENT VARIABLES ===================================================================== \n");
	_print(list.element);
	eslog::info(" = PER ELEMENT NODES VARIABLES =============================================================== \n");
	_print(list.enodes);
	eslog::info(" = PER ELEMENT GPS VARIABLES ================================================================= \n");
	_print(list.egps);
	eslog::info(" = PER NODE VARIABLES ======================================================================== \n");
	_print(list.node);
	eslog::info(" ============================================================================================= \n");
}

void Variable::clear()
{
	auto remove = [] (std::map<std::string, Variable*> &map) {
		for (auto it = map.begin(); it != map.end(); ++it) {
			if (it->second) {
				delete it->second;
			}
		}
	};

	remove(list.global);
	remove(list.element);
	remove(list.enodes);
	remove(list.egps);
	remove(list.node);
}

RangeVariable::RangeVariable(ECFRange &range)
: previous(range.value == 0 ? 1 : 0), range(range)
{

}

void RangeVariable::set(size_t interval, Evaluator::Params::General &param) const
{
	param.increment = 0;
	param.offset = 0;
	param.val = &range.value;
}

int RangeVariable::update(size_t interval) const
{
	return previous != range.value;
}

int RangeVariable::isconst(size_t interval) const
{
	return true;
}

void RangeVariable::updated()
{
	previous = range.value;
}

TimeVariable::TimeVariable(step::Time &time)
: previous(time.current == 0 ? 1 : 0), time(time)
{

}

void TimeVariable::set(size_t interval, Evaluator::Params::General &param) const
{
	param.increment = 0;
	param.offset = 0;
	param.val = &time.current;
}

int TimeVariable::update(size_t interval) const
{
	return previous != time.current;
}

int TimeVariable::isconst(size_t interval) const
{
	return true;
}

void TimeVariable::updated()
{
	previous = time.current;
}

FrequencyVariable::FrequencyVariable(step::Frequency &frequency)
: previous(frequency.current == 0 ? 1 : 0), frequency(frequency)
{

}

void FrequencyVariable::set(size_t interval, Evaluator::Params::General &param) const
{
	param.increment = 0;
	param.offset = 0;
	param.val = &frequency.current;
}

int FrequencyVariable::update(size_t interval) const
{
	return previous != frequency.current;
}

int FrequencyVariable::isconst(size_t interval) const
{
	return true;
}

void FrequencyVariable::updated()
{
	previous = frequency.current;
}

OutputVariable::OutputVariable(NamedData *data, int offset, int size)
: data(data), offset(offset), size(size)
{

}

void OutputVariable::set(size_t interval, Evaluator::Params::General &param) const
{
	param.increment = size;
	param.offset = offset;
	param.val = data->data.data();
}

int OutputVariable::update(size_t interval) const
{
	return data->updated;
}

int OutputVariable::isconst(size_t interval) const
{
	return false;
}

void OutputVariable::updated()
{
	data->updated = false;
}


SerializedEdataVariable::SerializedEdataVariable(serializededata<esint, double> *data, int offset, int size)
: data(data), offset(offset), size(size)
{

}

void SerializedEdataVariable::set(size_t interval, Evaluator::Params::General &param) const
{
	param.increment = size;
	param.offset = offset;
	param.val = data->begin(interval)->data();
}

int SerializedEdataVariable::update(size_t interval) const
{
	return false; // implement updating mechanism
}

int SerializedEdataVariable::isconst(size_t interval) const
{
	return false;
}

void SerializedEdataVariable::updated()
{

}

SerializedPointsVariable::SerializedPointsVariable(serializededata<esint, Point> *data, int offset)
: data(data), offset(offset)
{

}

void SerializedPointsVariable::set(size_t interval, Evaluator::Params::General &param) const
{
	param.increment = 3;
	param.offset = offset;
	param.val = &data->begin(interval)->data()->x;
}

int SerializedPointsVariable::update(size_t interval) const
{
	return false;
}

int SerializedPointsVariable::isconst(size_t interval) const
{
	return false;
}

void SerializedPointsVariable::updated()
{

}

ParameterVariable::ParameterVariable(serializededata<esint, double> *data, std::vector<int> &isconst, std::vector<int> &update, int offset, int size)
: data(data), constness(isconst), updating(update), offset(offset), size(size)
{

}

void ParameterVariable::set(size_t interval, Evaluator::Params::General &param) const
{
	param.increment = size;
	param.offset = offset;
	param.val = (data->begin() + interval)->data();
}

int ParameterVariable::update(size_t interval) const
{
	return updating[interval];
}

int ParameterVariable::isconst(size_t interval) const
{
	return constness[interval];
}

void ParameterVariable::updated()
{

}
