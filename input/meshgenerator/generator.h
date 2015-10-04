
#ifndef INPUT_MESHGENERATOR_GENERATOR_H_
#define INPUT_MESHGENERATOR_GENERATOR_H_

#include "elements/elements.h"

namespace esinput {

class Generator {

public:

	~Generator() { };

protected:

	Generator(): _processes(1) { }

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

};

class MeshGenerator {
public:
	MeshGenerator(int argc, char** argv);

	~MeshGenerator()
	{
		delete _generator;
	}

private:

	Generator *_generator;
};

}


#endif /* INPUT_MESHGENERATOR_GENERATOR_H_ */
