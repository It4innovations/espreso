
#include "evaluator.h"
#include "constevaluator.h"
#include "expressionevaluator.h"
#include "basis/utilities/parser.h"
#include "esinfo/envinfo.h"

#include "wrappers/exprtk/exprtk.h"

using namespace espreso;

Evaluator* Evaluator::create(const std::string &expression)
{
    if (expression == "") {
        return new Evaluator();
    }
//    if (StringCompare::contains(value, { "TABULAR" })) {
//        std::string value = Parser::strip(value.substr(value.find_first_of("[")));
//        value = value.substr(1, value.size() - 3);
//        std::vector<std::string> lines = Parser::split(value, ";");
//        std::vector<std::pair<double, double> > table;
//
//        for (size_t i = 0; i < lines.size(); i++) {
//            if (lines[i].size() == 0) {
//                continue;
//            }
//            std::vector<std::string> line = Parser::split(lines[i], ",");
//            if (line.size() != 2) {
//                eslog::globalerror("Invalid TABULAR data: %s\n", value.c_str());
//            }
//            table.push_back(std::make_pair(std::stod(line[0]), std::stod(line[1])));
//        }
//        evaluator = new TableInterpolationEvaluator(table);
//        evaluator->isset = isset;
//        return;
//    }
    std::vector<EvaluatorParameter> parameters;
    Exprtk::check(expression, parameters);
    if (parameters.size()) {
        return new ExpressionEvaluator(expression);
    } else {
        return new ConstEvaluator(Exprtk(expression, parameters).evaluate());
    }
}

Evaluator::Evaluator()
: parameters(info::env::threads, { EvaluatorParameter("dummy") })
{

}

double& Evaluator::getParameter(const std::string &name, int t)
{
    for (size_t i = 0; i < parameters[t].size(); ++i) {
        if (StringCompare::caseInsensitiveEq(parameters[t][i].name, name)) {
            return parameters[t][i].value;
        }
    }
    return parameters[t][0].value; // return dummy value
}
