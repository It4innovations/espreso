
#include "expressionevaluator.h"
#include "esinfo/envinfo.h"

using namespace espreso;

ExpressionEvaluator::ExpressionEvaluator(const std::string &expression)
: _expr(expression.c_str())
{
    _expression.resize(info::env::threads);
    #pragma omp parallel for
    for (int t = 0; t < info::env::threads; t++) {
        _expression[t] = new Exprtk(expression, parameters[t]);
    }
}

ExpressionEvaluator::~ExpressionEvaluator()
{
    for (int t = 0; t < info::env::threads; t++) {
        delete _expression[t];
    }
}




