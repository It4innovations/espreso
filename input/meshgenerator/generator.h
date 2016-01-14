
#ifndef INPUT_MESHGENERATOR_GENERATOR_H_
#define INPUT_MESHGENERATOR_GENERATOR_H_

#include "../loader.h"
#include "elements/elements.h"
#include "settings.h"

namespace esinput {

class Generator: public InternalLoader {

public:
	Generator(int argc, char** argv, size_t index, size_t size): _settings(argc, argv, index, size) { };
	Generator(const Settings &settings): _settings(settings) { };

	virtual void points(mesh::Coordinates &coordinates) = 0;
	virtual void elements(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts) = 0;
	virtual void fixPoints(std::vector<std::vector<eslocal> > &fixPoints) = 0;
	virtual void boundaryConditions(mesh::Coordinates &coordinates) = 0;
	virtual void corners(mesh::Boundaries &boundaries) = 0;
	virtual void clusterBoundaries(mesh::Boundaries &boundaries) = 0;

	bool manualPartition()
	{
		return _settings.useMetis;
	}

	virtual ~Generator() { };

protected:
	const Settings _settings;
};

class MeshGenerator: public InternalLoader {
public:
	MeshGenerator(int argc, char** argv, size_t index, size_t size);
	MeshGenerator(Generator *generator): _generator(generator) { };

	bool manualPartition();

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts);
	void fixPoints(std::vector<std::vector<eslocal> > &fixPoints);
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
