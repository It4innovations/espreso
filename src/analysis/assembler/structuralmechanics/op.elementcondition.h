
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ELEMENTCONDITION_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ELEMENTCONDITION_H_

#include "analysis/assembler/general/boundarycondition.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct ElementCondition: BoundaryCondition {

	ElementCondition()
	: behaviour(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN)
	{

	}

	void activate(ECFExpressionVector *expressionVector, StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
	{
		this->behaviour = behaviour;
		BoundaryCondition::activate(expressionVector);
	}

	void activate(ECFExpression *expression, StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
	{
		this->behaviour = behaviour;
		BoundaryCondition::activate(expression);
	}

	StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ELEMENTCONDITION_H_ */
