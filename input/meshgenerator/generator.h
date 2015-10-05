
#ifndef INPUT_MESHGENERATOR_GENERATOR_H_
#define INPUT_MESHGENERATOR_GENERATOR_H_

#include "../loader.h"
#include "elements/elements.h"

namespace esinput {

class Generator {

public:

	virtual void points(mesh::Coordinates &data) = 0;
	virtual void elements(std::vector<mesh::Element*> &data) = 0;

	bool manualPostProcessing()
	{
		return false;
	}

	virtual ~Generator() { };

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

class MeshGenerator: public InternalLoader {
public:
	MeshGenerator(int argc, char** argv, int rank, int size);

	void points(mesh::Coordinates &data);
	void elements(std::vector<mesh::Element*> &data);

	~MeshGenerator()
	{
		delete _generator;
	}

private:

	Generator *_generator;
};

}


#endif /* INPUT_MESHGENERATOR_GENERATOR_H_ */
