
#ifndef SRC_BASIS_EVALUATOR_TABLEINTERPOLATIONEVALUATOR_H_
#define SRC_BASIS_EVALUATOR_TABLEINTERPOLATIONEVALUATOR_H_

#include "evaluator.h"

#include <utility>
#include <vector>

namespace espreso {

class TableInterpolationEvaluator: public Evaluator {

public:
	TableInterpolationEvaluator(const std::vector<std::pair<double, double> > &table);

	double evaluate() const;

protected:
	std::vector<std::pair<double, double> > _table;
};

}



#endif /* SRC_BASIS_EVALUATOR_TABLEINTERPOLATIONEVALUATOR_H_ */
