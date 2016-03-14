
#ifndef INPUT_MESHGENERATOR_GENERATOR_H_
#define INPUT_MESHGENERATOR_GENERATOR_H_

#include "../loader.h"
#include "elements/elements.h"
#include "settings.h"

namespace espreso {
namespace input {

class Generator: public InternalLoader {

protected:
	Generator(const Settings &settings): _settings(settings) { };
	virtual ~Generator() { };

	void elements(std::vector<Element*> &elements, std::vector<eslocal> &parts)
	{
		elementsMesh(elements, parts);
		elementsMaterials(elements, parts);
	}

	virtual void elementsMesh(std::vector<Element*> &elements, std::vector<eslocal> &parts) = 0;
	virtual void elementsMaterials(std::vector<Element*> &elements, std::vector<eslocal> &parts) = 0;

	virtual void points(Coordinates &coordinates) = 0;
	virtual void fixPoints(std::vector<std::vector<eslocal> > &fixPoints) = 0;
	virtual void boundaryConditions(Coordinates &coordinates) = 0;
	virtual void corners(Boundaries &boundaries) = 0;
	virtual void clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours) = 0;

	bool manualPartition()
	{
		return _settings.useMetis;
	}

	const Settings _settings;
};

}
}


#endif /* INPUT_MESHGENERATOR_GENERATOR_H_ */
