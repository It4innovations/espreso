
#ifndef SRC_BASIS_EVALUATOR_TABLEINTERPOLATIONEVALUATOR_H_
#define SRC_BASIS_EVALUATOR_TABLEINTERPOLATIONEVALUATOR_H_

#include "evaluator.h"

#include <utility>
#include <vector>

namespace espreso {

class TableInterpolationEvaluator: public Evaluator {

public:
	TableInterpolationEvaluator(const std::vector<std::pair<double, double> > &table);

	Type type() { return Type::TABLE_INTERPOLATION; }
	Evaluator* copy() const { return new TableInterpolationEvaluator(*this); }

	void evalVector(esint size, esint increment, const Params &params, double *results) const;
	void evalFiltered(esint size, esint increment, const esint *elements, const esint *distribution, const Params &params, double *results) const;
	void evalSelectedSparse(esint size, esint increment, const esint *selection, const Params &params, double *results) const;
	void evalSelectedDense(esint size, esint increment, const esint *selection, const Params &params, double *results) const;

	std::string getEXPRTKForm() const;

protected:
	std::vector<std::pair<double, double> > _table;
};

}



#endif /* SRC_BASIS_EVALUATOR_TABLEINTERPOLATIONEVALUATOR_H_ */
