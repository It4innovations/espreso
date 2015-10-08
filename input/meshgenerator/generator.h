
#ifndef INPUT_MESHGENERATOR_GENERATOR_H_
#define INPUT_MESHGENERATOR_GENERATOR_H_

#include "../loader.h"
#include "elements/elements.h"

namespace esinput {

class Generator {

public:
	virtual void points(mesh::Coordinates &coordinates) = 0;
	virtual void elements(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts) = 0;
	virtual void fixPoints(std::vector<eslocal> &fixPoints) = 0;
	virtual void boundaryConditions(mesh::Coordinates &coordinates) = 0;
	virtual void corners(mesh::Boundaries &boundaries) = 0;
	virtual void clusterBoundaries(mesh::Boundaries &boundaries) = 0;

	virtual ~Generator() { };
};

class MeshGenerator: public InternalLoader {
public:
	MeshGenerator(int argc, char** argv, int rank, int size);

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts);
	void fixPoints(std::vector<eslocal> &fixPoints);
	void boundaryConditions(mesh::Coordinates &coordinates);
	void corners(mesh::Boundaries &boundaries);
	void clusterBoundaries(mesh::Boundaries &boundaries);

	~MeshGenerator()
	{
		delete _generator;
	}

private:

	Generator *_generator;
};

}


#endif /* INPUT_MESHGENERATOR_GENERATOR_H_ */
