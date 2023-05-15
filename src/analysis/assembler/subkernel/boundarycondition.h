
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BOUNDARYCONDITION_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BOUNDARYCONDITION_H_

#include "subkernel.h"

namespace espreso {

struct BoundaryCondition: public SubKernel {
	ECFExpression *expression;
	ECFExpressionVector *expressionVector;
	double *rhs;

	BoundaryCondition()
	: expression(nullptr), expressionVector(nullptr), rhs(nullptr)
	{
		isconst = false;
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE | Assembler::ITERATION;
	}

	void activate(ECFExpressionVector *expressionVector, double *rhs)
	{
		this->expressionVector = expressionVector;
		this->rhs = rhs;
		if (this->expressionVector) {
			this->isactive = 1;
		}
	}

	void activate(ECFExpression *expression, double *rhs)
	{
		this->expression = expression;
		this->rhs = rhs;
		if (this->expression) {
			this->isactive = 1;
		}
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BOUNDARYCONDITION_H_ */
