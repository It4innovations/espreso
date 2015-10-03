
#ifndef INPUT_MESHGENERATOR_MESHGENERATOR_H_
#define INPUT_MESHGENERATOR_MESHGENERATOR_H_

#include "elements/elements.h"

namespace esinput {

class MeshGenerator {

public:

	MeshGenerator(const std::string &configFile): _processes(1)
	{
		_generator = NULL;
	}

	~MeshGenerator() { }

protected:
	size_t assumedProcessCount()
	{
		return _processes;
	}

	void fillCluster(int rank, size_t cluster[]);

	void mesh(mesh::Mesh &mesh, const size_t cluster[]);

	void setDirichlet(mesh::Mesh &mesh, const size_t cluster[], size_t dirichlet);
	void setForces(mesh::Mesh &mesh, const size_t cluster[]);

	void fillGlobalBoundaries(mesh::Boundaries &boundaries, const size_t cluster[]);

	void setFixPoints(mesh::Mesh &mesh);
	void setCorners(
			mesh::Boundaries &boundaries,
			const size_t number[],
			const bool corners,
			const bool edges,
			const bool surface);

	size_t _processes;
	MeshGenerator *_generator;
};

}


#endif /* INPUT_MESHGENERATOR_MESHGENERATOR_H_ */
