
#ifndef INPUT_MESHGENERATOR_GENERATOR_H_
#define INPUT_MESHGENERATOR_GENERATOR_H_

#include "../loader.h"
#include "elements/elements.h"
#include "settings.h"

namespace esinput {

class Generator: public InternalLoader {

protected:
	Generator(const Settings &settings): _settings(settings) { };
	virtual ~Generator() { };

	void elements(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts)
	{
		elementsMesh(elements, parts);
		elementsMaterials(elements, parts);
	}

	virtual void elementsMesh(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts) = 0;
	virtual void elementsMaterials(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts) = 0;

	virtual void points(mesh::Coordinates &coordinates) = 0;
	virtual void fixPoints(std::vector<std::vector<eslocal> > &fixPoints) = 0;
	virtual void boundaryConditions(mesh::Coordinates &coordinates) = 0;
	virtual void corners(mesh::Boundaries &boundaries) = 0;
	virtual void clusterBoundaries(mesh::Boundaries &boundaries) = 0;

	bool manualPartition()
	{
		return _settings.useMetis;
	}

	const Settings _settings;
};

}


#endif /* INPUT_MESHGENERATOR_GENERATOR_H_ */
