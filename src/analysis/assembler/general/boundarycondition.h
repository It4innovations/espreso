
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BOUNDARYCONDITION_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BOUNDARYCONDITION_H_

#include "subkernel.h"
#include "basis/evaluator/evaluator.h"
#include "config/holders/expression.h"

namespace espreso {

struct BoundaryCondition: public SubKernel {
	ECFExpression *expression;
	ECFExpressionVector *expressionVector;

	BoundaryCondition()
	: expression(nullptr), expressionVector(nullptr)
	{
		isconst = false;
		action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION;
	}

	void activate(ECFExpressionVector *expressionVector)
	{
		this->expressionVector = expressionVector;
		if (this->expressionVector) {
			this->isconst = expressionVector->x.evaluator->isConst() && expressionVector->y.evaluator->isConst() && expressionVector->z.evaluator->isConst();
			this->isactive = 1;
			this->needCoordinates = expressionVector->x.evaluator->needCoordinates() && expressionVector->y.evaluator->needCoordinates() && expressionVector->z.evaluator->needCoordinates();
		}
	}

	void activate(ECFExpression *expression)
	{
		this->expression = expression;
		if (this->expression) {
			this->isconst = expression->evaluator->isConst();
			this->isactive = 1;
			this->needCoordinates = expression->evaluator->needCoordinates();
		}
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BOUNDARYCONDITION_H_ */
