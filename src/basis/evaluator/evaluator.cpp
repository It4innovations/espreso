
#include "evaluator.h"
#include "constevaluator.h"
#include "expressionevaluator.h"

#include "wrappers/exprtk/exprtk.h"

using namespace espreso;

Evaluator* Evaluator::create(const std::string &expression)
{
	if (expression == "") {
		return new Evaluator();
	}
//	if (StringCompare::contains(value, { "TABULAR" })) {
//		std::string value = Parser::strip(value.substr(value.find_first_of("[")));
//		value = value.substr(1, value.size() - 3);
//		std::vector<std::string> lines = Parser::split(value, ";");
//		std::vector<std::pair<double, double> > table;
//
//		for (size_t i = 0; i < lines.size(); i++) {
//			if (lines[i].size() == 0) {
//				continue;
//			}
//			std::vector<std::string> line = Parser::split(lines[i], ",");
//			if (line.size() != 2) {
//				eslog::globalerror("Invalid TABULAR data: %s\n", value.c_str());
//			}
//			table.push_back(std::make_pair(std::stod(line[0]), std::stod(line[1])));
//		}
//		evaluator = new TableInterpolationEvaluator(table);
//		evaluator->isset = isset;
//		return;
//	}
	std::vector<EvaluatorParameter> parameters;
	Exprtk::check(expression, parameters);
	if (parameters.size()) {
		return new ExpressionEvaluator(expression);
	} else {
		return new ConstEvaluator(Exprtk(expression, parameters).evaluate());
	}
}

