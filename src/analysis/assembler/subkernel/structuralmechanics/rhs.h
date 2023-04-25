
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_RHS_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_RHS_H_

#include "subkernel.h"

namespace espreso {

struct StructuralMechanicsRHS: public SubKernel {
	ECFExpression *expression;
	ECFExpressionVector *expressionVector;
	double *rhs;

	StructuralMechanicsRHS()
	: expression(nullptr), expressionVector(nullptr), rhs(nullptr)
	{
		isconst = false;
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE;
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

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_RHS_H_ */
