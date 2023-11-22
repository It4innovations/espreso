
#include "analysis/assembler/kernel/structuralmechanics/kernel.h"
#include "analysis/assembler/module/structuralmechanics.h"

namespace espreso {

template <>
void StructuralMechanics::runGroup<Element::CODE::TETRA4>(Action action, size_t interval, StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
{
	runElasticity<Element::CODE::TETRA4, 4, StructuralMechanicsGPC::TETRA4, 3, 3, Behaviour::VOLUME>(subkernels[interval], action);
}

}
