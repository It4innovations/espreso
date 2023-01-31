
#ifndef SRC_BASIS_EVALUATOR_TABLEEVALUATOR_H_
#define SRC_BASIS_EVALUATOR_TABLEEVALUATOR_H_

#include "evaluator.h"

#include <cstddef>
#include <vector>

namespace espreso {

class TableEvaluator: public Evaluator {

public:
	enum class TableProperty {
		TEMPERATURE,
		TIME
	};

	TableEvaluator(
			const std::vector<std::vector<std::vector<double> > > &table,
			const std::vector<TableProperty> &properties,
			const std::vector<std::vector<double> > &axis);

	virtual Evaluator* copy() const { return new TableEvaluator(*this); }

	void evalVectorInc(esint size, esint increment, const Params &params, double *results) const;
	void evalFilteredInc(esint size, esint increment, const esint *elements, const esint *distribution, const Params &params, double *results) const;
	void evalSelectedSparseInc(esint size, esint increment, const esint *selection, const Params &params, double *results) const;
	void evalSelectedDenseInc(esint size, esint increment, const esint *selection, const Params &params, double *results) const;

	std::string getEXPRTKForm() const;
	std::string toString() const { return "TABLE"; }

protected:
	size_t _dimension;
	std::vector<std::vector<std::vector<double> > > _table;
	std::vector<TableProperty> _properties;
	std::vector<std::vector<double> > _axis;
};

}

#endif /* SRC_BASIS_EVALUATOR_TABLEEVALUATOR_H_ */
