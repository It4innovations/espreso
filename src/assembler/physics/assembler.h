
#ifndef SRC_ASSEMBLER_PHYSICS_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_ASSEMBLER_H_

#include "esmesh.h"
#include "essolver.h"

namespace espreso {

struct Physics {

	virtual void assemble() =0;
	virtual void save() =0;

	size_t DOFs;

	Physics(const Mesh &mesh, size_t DOFs): _mesh(mesh), DOFs(DOFs) {};
	virtual ~Physics() {};

protected:
	const Mesh& _mesh;
};

}

#endif /* SRC_ASSEMBLER_PHYSICS_ASSEMBLER_H_ */
