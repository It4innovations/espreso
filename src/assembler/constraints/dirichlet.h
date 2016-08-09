
#ifndef SRC_ASSEMBLER_CONSTRAINTS_DIRICHLET_H_
#define SRC_ASSEMBLER_CONSTRAINTS_DIRICHLET_H_

#include "constraints.h"

namespace espreso {

class Dirichlet: public ConstraintsBase
{
public:
	Dirichlet(Mesh &mesh, Physics &physics): ConstraintsBase(mesh, physics) {};

	virtual void insertDirichletToB1(const std::vector<Element*> &nodes, const Coordinates &coordinates, const std::vector<Property> &DOFs);

};

}

#endif /* SRC_ASSEMBLER_CONSTRAINTS_DIRICHLET_H_ */
