
#ifndef SRC_ASSEMBLER_CONSTRAINTS_INEQUALITYCONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_INEQUALITYCONSTRAINTS_H_

#include "constraints.h"

namespace espreso {

struct InequalityConstraints
{
	static void insertLowerBoundToB1(Constraints &constraints, const std::vector<Property> &eDOFs, const std::vector<Property> &boundDOFs);
};

}



#endif /* SRC_ASSEMBLER_CONSTRAINTS_INEQUALITYCONSTRAINTS_H_ */
