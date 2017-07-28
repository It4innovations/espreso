
#ifndef SRC_ASSEMBLER_PHYSICS_BOUNDARYBASED3D_H_
#define SRC_ASSEMBLER_PHYSICS_BOUNDARYBASED3D_H_

#include "physics.h"

namespace espreso {

struct BoundaryBased3D: virtual public Physics
{
	BoundaryBased3D();

	void extractBoundaryNodes();
	void boundaryTriangularization(std::vector<eslocal> &elements, std::vector<double> &coordinates, size_t domain);

protected:
	std::vector<std::vector<eslocal> > _boundaryIndices;
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_BOUNDARYBASED3D_H_ */
