
#include "analysis/assembler/kernel/structuralmechanics/kernel.h"
#include "analysis/assembler/module/structuralmechanics.h"

namespace espreso {

template <>
void StructuralMechanics::runGroup<Element::CODE::SQUARE4>(Action action, size_t interval, StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
{
	if (behaviour == StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::AXISYMMETRIC) {
		runElasticity<Element::CODE::SQUARE4, 4, StructuralMechanicsGPC::SQUARE4, 2, 2, Behaviour::AXISYMMETRIC>(subkernels[interval], action);
	} else {
		runElasticity<Element::CODE::SQUARE4, 4, StructuralMechanicsGPC::SQUARE4, 2, 2, Behaviour::PLANE>(subkernels[interval], action);
	}
}

}
