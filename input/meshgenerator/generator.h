
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

	virtual ~Generator() { };

protected:

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

};

class MeshGenerator: public InternalLoader {
public:
	MeshGenerator(int argc, char** argv, int rank, int size);

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts);
	void fixPoints(std::vector<eslocal> &fixPoints);

	~MeshGenerator()
	{
		delete _generator;
	}

private:

	Generator *_generator;
};

}


#endif /* INPUT_MESHGENERATOR_GENERATOR_H_ */
