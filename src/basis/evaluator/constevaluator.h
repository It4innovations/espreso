
#ifndef SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_
#define SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_

#include "evaluator.h"

namespace espreso {

class ConstEvaluator: public Evaluator {

public:
    ConstEvaluator(double value): _expr(std::to_string(value)), _value(value) {}

    double evaluate(int t) const { return _value; }

    const char* expression() const { return _expr.c_str(); }
protected:
    std::string _expr;
    const double _value;
};

}


#endif /* SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_ */
