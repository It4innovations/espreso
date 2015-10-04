
#ifndef GENERATOR_H_
#define GENERATOR_H_

#include "esmesh.h"
#include "../settings.h"
#include "../elements/elements.h"
#include "../utils.h"

namespace esinput {

class Generator {

public:
	eslocal assumedProcessCount()
	{
		return _processes;
	}

	virtual void fillCluster(int rank, size_t cluster[]) = 0;

	virtual void mesh(mesh::Mesh &mesh, const size_t cluster[]) = 0;

	virtual void setDirichlet(mesh::Mesh &mesh, const size_t cluster[], size_t dirichlet) = 0;
	virtual void setForces(mesh::Mesh &mesh, const size_t cluster[]) = 0;

	virtual void fillGlobalBoundaries(mesh::Boundaries &boundaries, const size_t cluster[]) = 0;

	virtual void setFixPoints(mesh::Mesh &mesh) = 0;
	virtual void setCorners(
			mesh::Boundaries &boundaries,
			const size_t number[],
			const bool corners,
			const bool edges,
			const bool surface) = 0;

	virtual ~Generator() { }

protected:
	MeshGenerator(): _processes(1) { };

	eslocal _processes;
};

}


#endif /* GENERATOR_H_ */
