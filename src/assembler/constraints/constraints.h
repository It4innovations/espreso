
#ifndef SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_

#include "esmesh.h"
#include "essolver.h"

namespace espreso {

class ConstraintsBase
{
public:
	ConstraintsBase(Mesh &mesh, Physics &physics): _mesh(mesh), _physics(physics) {};

	Mesh &_mesh;
	Physics &_physics;
};

}


#endif /* SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_ */
